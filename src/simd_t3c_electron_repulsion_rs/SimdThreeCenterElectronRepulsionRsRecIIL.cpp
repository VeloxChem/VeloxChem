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


#include "SimdThreeCenterElectronRepulsionRsRecIIL.hpp"

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
compute_rs_iil_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_iil_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 1154910, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 5746 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 1154910, 837208, 47929, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 20,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 28, 3, 20,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 101, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 104, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 107, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 110, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 116, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 119, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 125, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 131, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 134, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 137, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 140, 0, 3, 39, 40,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 40, 41,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 146, 0, 3, 41, 42,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 149, 0, 3, 42, 43,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 43, 44,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 155, 0, 3, 44, 45,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 45, 46,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 161, 0, 3, 46, 47,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 164, 0, 3, 47, 48,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 167, 0, 3, 48, 49,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 170, 0, 3, 7, 8,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 176, 0, 3, 8, 9,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 9, 10,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 10, 11,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 194, 0, 3, 11, 12,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 200, 0, 3, 12, 13,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 206, 0, 3, 13, 14,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 14, 15,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 15, 16,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 224, 0, 3, 16, 17,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 230, 0, 3, 17, 18,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 236, 0, 3, 18, 19,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 19, 20,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 20, 21,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 254, 0, 3, 21, 22,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 260, 0, 3, 22, 23,
                                                                       95, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 266, 0, 3, 23, 24,
                                                                       98, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 24, 25,
                                                                       101, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 25, 26,
                                                                       104, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 284, 0, 3, 29, 30,
                                                                       110, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 290, 0, 3, 30, 31,
                                                                       113, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 296, 0, 3, 31, 32,
                                                                       116, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 302, 0, 3, 32, 33,
                                                                       119, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 33, 34,
                                                                       122, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 314, 0, 3, 34, 35,
                                                                       125, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 320, 0, 3, 35, 36,
                                                                       128, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 326, 0, 3, 36, 37,
                                                                       131, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 332, 0, 3, 37, 38,
                                                                       134, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 38, 39,
                                                                       137, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 344, 0, 3, 39, 40,
                                                                       140, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 350, 0, 3, 40, 41,
                                                                       143, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 356, 0, 3, 41, 42,
                                                                       146, 149, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 362, 0, 3, 42, 43,
                                                                       149, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 43, 44,
                                                                       152, 155, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 374, 0, 3, 44, 45,
                                                                       155, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 380, 0, 3, 45, 46,
                                                                       158, 161, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 386, 0, 3, 46, 47,
                                                                       161, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 392, 0, 3, 47, 48,
                                                                       164, 167, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 50, 53,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 53, 56,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 56, 59,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 59, 62,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 62, 65,
                                                                       194, 200, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 65, 68,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 68, 71,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 71, 74,
                                                                       212, 218, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 74, 77,
                                                                       218, 224, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 77, 80,
                                                                       224, 230, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 80, 83,
                                                                       230, 236, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 83, 86,
                                                                       236, 242, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 86, 89,
                                                                       242, 248, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 528, 0, 3, 89, 92,
                                                                       248, 254, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 538, 0, 3, 92, 95,
                                                                       254, 260, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 548, 0, 3, 95, 98,
                                                                       260, 266, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 558, 0, 3, 98,
                                                                       101, 266, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 568, 0, 3, 101,
                                                                       104, 272, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 578, 0, 3, 110,
                                                                       113, 284, 290, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 588, 0, 3, 113,
                                                                       116, 290, 296, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 598, 0, 3, 116,
                                                                       119, 296, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 608, 0, 3, 119,
                                                                       122, 302, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 618, 0, 3, 122,
                                                                       125, 308, 314, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 628, 0, 3, 125,
                                                                       128, 314, 320, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 128,
                                                                       131, 320, 326, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 648, 0, 3, 131,
                                                                       134, 326, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 658, 0, 3, 134,
                                                                       137, 332, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 668, 0, 3, 137,
                                                                       140, 338, 344, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 678, 0, 3, 140,
                                                                       143, 344, 350, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 688, 0, 3, 143,
                                                                       146, 350, 356, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 698, 0, 3, 146,
                                                                       149, 356, 362, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 708, 0, 3, 149,
                                                                       152, 362, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 718, 0, 3, 152,
                                                                       155, 368, 374, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 728, 0, 3, 155,
                                                                       158, 374, 380, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 738, 0, 3, 158,
                                                                       161, 380, 386, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 748, 0, 3, 161,
                                                                       164, 386, 392, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 758, 0, 3, 170,
                                                                       176, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 773, 0, 3, 176,
                                                                       182, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 788, 0, 3, 182,
                                                                       188, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 803, 0, 3, 188,
                                                                       194, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 818, 0, 3, 194,
                                                                       200, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 833, 0, 3, 200,
                                                                       206, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 848, 0, 3, 206,
                                                                       212, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 863, 0, 3, 212,
                                                                       218, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 878, 0, 3, 218,
                                                                       224, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 893, 0, 3, 224,
                                                                       230, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 908, 0, 3, 230,
                                                                       236, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 923, 0, 3, 236,
                                                                       242, 508, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 938, 0, 3, 242,
                                                                       248, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 953, 0, 3, 248,
                                                                       254, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 968, 0, 3, 254,
                                                                       260, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 983, 0, 3, 260,
                                                                       266, 548, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 998, 0, 3, 266,
                                                                       272, 558, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1013, 0, 3, 284,
                                                                       290, 578, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1028, 0, 3, 290,
                                                                       296, 588, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1043, 0, 3, 296,
                                                                       302, 598, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 302,
                                                                       308, 608, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1073, 0, 3, 308,
                                                                       314, 618, 628, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1088, 0, 3, 314,
                                                                       320, 628, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1103, 0, 3, 320,
                                                                       326, 638, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1118, 0, 3, 326,
                                                                       332, 648, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1133, 0, 3, 332,
                                                                       338, 658, 668, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1148, 0, 3, 338,
                                                                       344, 668, 678, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1163, 0, 3, 344,
                                                                       350, 678, 688, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1178, 0, 3, 350,
                                                                       356, 688, 698, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1193, 0, 3, 356,
                                                                       362, 698, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1208, 0, 3, 362,
                                                                       368, 708, 718, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1223, 0, 3, 368,
                                                                       374, 718, 728, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1238, 0, 3, 374,
                                                                       380, 728, 738, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1253, 0, 3, 380,
                                                                       386, 738, 748, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 398,
                                                                       408, 758, 773, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1289, 0, 3, 408,
                                                                       418, 773, 788, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1310, 0, 3, 418,
                                                                       428, 788, 803, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1331, 0, 3, 428,
                                                                       438, 803, 818, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 438,
                                                                       448, 818, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1373, 0, 3, 448,
                                                                       458, 833, 848, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1394, 0, 3, 458,
                                                                       468, 848, 863, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1415, 0, 3, 468,
                                                                       478, 863, 878, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 478,
                                                                       488, 878, 893, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1457, 0, 3, 488,
                                                                       498, 893, 908, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1478, 0, 3, 498,
                                                                       508, 908, 923, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1499, 0, 3, 508,
                                                                       518, 923, 938, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 518,
                                                                       528, 938, 953, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1541, 0, 3, 528,
                                                                       538, 953, 968, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1562, 0, 3, 538,
                                                                       548, 968, 983, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1583, 0, 3, 548,
                                                                       558, 983, 998, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 578,
                                                                       588, 1013, 1028, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1625, 0, 3, 588,
                                                                       598, 1028, 1043, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1646, 0, 3, 598,
                                                                       608, 1043, 1058, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1667, 0, 3, 608,
                                                                       618, 1058, 1073, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 618,
                                                                       628, 1073, 1088, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1709, 0, 3, 628,
                                                                       638, 1088, 1103, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1730, 0, 3, 638,
                                                                       648, 1103, 1118, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1751, 0, 3, 648,
                                                                       658, 1118, 1133, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 658,
                                                                       668, 1133, 1148, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1793, 0, 3, 668,
                                                                       678, 1148, 1163, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1814, 0, 3, 678,
                                                                       688, 1163, 1178, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1835, 0, 3, 688,
                                                                       698, 1178, 1193, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 698,
                                                                       708, 1193, 1208, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1877, 0, 3, 708,
                                                                       718, 1208, 1223, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1898, 0, 3, 718,
                                                                       728, 1223, 1238, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1919, 0, 3, 728,
                                                                       738, 1238, 1253, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 758,
                                                                       773, 1268, 1289, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1968, 0, 3, 773,
                                                                       788, 1289, 1310, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1996, 0, 3, 788,
                                                                       803, 1310, 1331, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 803,
                                                                       818, 1331, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2052, 0, 3, 818,
                                                                       833, 1352, 1373, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 833,
                                                                       848, 1373, 1394, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 848,
                                                                       863, 1394, 1415, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2136, 0, 3, 863,
                                                                       878, 1415, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2164, 0, 3, 878,
                                                                       893, 1436, 1457, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2192, 0, 3, 893,
                                                                       908, 1457, 1478, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2220, 0, 3, 908,
                                                                       923, 1478, 1499, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2248, 0, 3, 923,
                                                                       938, 1499, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 938,
                                                                       953, 1520, 1541, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2304, 0, 3, 953,
                                                                       968, 1541, 1562, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2332, 0, 3, 968,
                                                                       983, 1562, 1583, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2360, 0, 3, 1013,
                                                                       1028, 1604, 1625, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2388, 0, 3, 1028,
                                                                       1043, 1625, 1646, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2416, 0, 3, 1043,
                                                                       1058, 1646, 1667, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2444, 0, 3, 1058,
                                                                       1073, 1667, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2472, 0, 3, 1073,
                                                                       1088, 1688, 1709, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2500, 0, 3, 1088,
                                                                       1103, 1709, 1730, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2528, 0, 3, 1103,
                                                                       1118, 1730, 1751, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2556, 0, 3, 1118,
                                                                       1133, 1751, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2584, 0, 3, 1133,
                                                                       1148, 1772, 1793, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2612, 0, 3, 1148,
                                                                       1163, 1793, 1814, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2640, 0, 3, 1163,
                                                                       1178, 1814, 1835, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2668, 0, 3, 1178,
                                                                       1193, 1835, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2696, 0, 3, 1193,
                                                                       1208, 1856, 1877, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2724, 0, 3, 1208,
                                                                       1223, 1877, 1898, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2752, 0, 3, 1223,
                                                                       1238, 1898, 1919, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2780, 0, 3, 1268,
                                                                       1289, 1940, 1968, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2816, 0, 3, 1289,
                                                                       1310, 1968, 1996, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2852, 0, 3, 1310,
                                                                       1331, 1996, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2888, 0, 3, 1331,
                                                                       1352, 2024, 2052, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2924, 0, 3, 1352,
                                                                       1373, 2052, 2080, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2960, 0, 3, 1373,
                                                                       1394, 2080, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2996, 0, 3, 1394,
                                                                       1415, 2108, 2136, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3032, 0, 3, 1415,
                                                                       1436, 2136, 2164, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3068, 0, 3, 1436,
                                                                       1457, 2164, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3104, 0, 3, 1457,
                                                                       1478, 2192, 2220, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3140, 0, 3, 1478,
                                                                       1499, 2220, 2248, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3176, 0, 3, 1499,
                                                                       1520, 2248, 2276, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3212, 0, 3, 1520,
                                                                       1541, 2276, 2304, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3248, 0, 3, 1541,
                                                                       1562, 2304, 2332, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3284, 0, 3, 1604,
                                                                       1625, 2360, 2388, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3320, 0, 3, 1625,
                                                                       1646, 2388, 2416, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3356, 0, 3, 1646,
                                                                       1667, 2416, 2444, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3392, 0, 3, 1667,
                                                                       1688, 2444, 2472, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3428, 0, 3, 1688,
                                                                       1709, 2472, 2500, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3464, 0, 3, 1709,
                                                                       1730, 2500, 2528, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3500, 0, 3, 1730,
                                                                       1751, 2528, 2556, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3536, 0, 3, 1751,
                                                                       1772, 2556, 2584, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3572, 0, 3, 1772,
                                                                       1793, 2584, 2612, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3608, 0, 3, 1793,
                                                                       1814, 2612, 2640, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3644, 0, 3, 1814,
                                                                       1835, 2640, 2668, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3680, 0, 3, 1835,
                                                                       1856, 2668, 2696, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3716, 0, 3, 1856,
                                                                       1877, 2696, 2724, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3752, 0, 3, 1877,
                                                                       1898, 2724, 2752, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3788, 0, 3, 1940,
                                                                       1968, 2780, 2816, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3833, 0, 3, 1968,
                                                                       1996, 2816, 2852, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3878, 0, 3, 1996,
                                                                       2024, 2852, 2888, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3923, 0, 3, 2024,
                                                                       2052, 2888, 2924, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 2052,
                                                                       2080, 2924, 2960, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4013, 0, 3, 2080,
                                                                       2108, 2960, 2996, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4058, 0, 3, 2108,
                                                                       2136, 2996, 3032, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4103, 0, 3, 2136,
                                                                       2164, 3032, 3068, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4148, 0, 3, 2164,
                                                                       2192, 3068, 3104, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4193, 0, 3, 2192,
                                                                       2220, 3104, 3140, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4238, 0, 3, 2220,
                                                                       2248, 3140, 3176, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4283, 0, 3, 2248,
                                                                       2276, 3176, 3212, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4328, 0, 3, 2276,
                                                                       2304, 3212, 3248, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4373, 0, 3, 2360,
                                                                       2388, 3284, 3320, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4418, 0, 3, 2388,
                                                                       2416, 3320, 3356, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4463, 0, 3, 2416,
                                                                       2444, 3356, 3392, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4508, 0, 3, 2444,
                                                                       2472, 3392, 3428, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4553, 0, 3, 2472,
                                                                       2500, 3428, 3464, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4598, 0, 3, 2500,
                                                                       2528, 3464, 3500, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4643, 0, 3, 2528,
                                                                       2556, 3500, 3536, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4688, 0, 3, 2556,
                                                                       2584, 3536, 3572, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4733, 0, 3, 2584,
                                                                       2612, 3572, 3608, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4778, 0, 3, 2612,
                                                                       2640, 3608, 3644, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4823, 0, 3, 2640,
                                                                       2668, 3644, 3680, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4868, 0, 3, 2668,
                                                                       2696, 3680, 3716, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4913, 0, 3, 2696,
                                                                       2724, 3716, 3752, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4958, 0, 3, 2780,
                                                                       2816, 3788, 3833, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5013, 0, 3, 2816,
                                                                       2852, 3833, 3878, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5068, 0, 3, 2852,
                                                                       2888, 3878, 3923, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5123, 0, 3, 2888,
                                                                       2924, 3923, 3968, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5178, 0, 3, 2924,
                                                                       2960, 3968, 4013, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5233, 0, 3, 2960,
                                                                       2996, 4013, 4058, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 2996,
                                                                       3032, 4058, 4103, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5343, 0, 3, 3032,
                                                                       3068, 4103, 4148, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5398, 0, 3, 3068,
                                                                       3104, 4148, 4193, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5453, 0, 3, 3104,
                                                                       3140, 4193, 4238, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5508, 0, 3, 3140,
                                                                       3176, 4238, 4283, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5563, 0, 3, 3176,
                                                                       3212, 4283, 4328, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5618, 0, 3, 3284,
                                                                       3320, 4373, 4418, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5673, 0, 3, 3320,
                                                                       3356, 4418, 4463, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5728, 0, 3, 3356,
                                                                       3392, 4463, 4508, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5783, 0, 3, 3392,
                                                                       3428, 4508, 4553, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5838, 0, 3, 3428,
                                                                       3464, 4553, 4598, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5893, 0, 3, 3464,
                                                                       3500, 4598, 4643, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5948, 0, 3, 3500,
                                                                       3536, 4643, 4688, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 6003, 0, 3, 3536,
                                                                       3572, 4688, 4733, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 6058, 0, 3, 3572,
                                                                       3608, 4733, 4778, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 6113, 0, 3, 3608,
                                                                       3644, 4778, 4823, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 6168, 0, 3, 3644,
                                                                       3680, 4823, 4868, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 6223, 0, 3, 3680,
                                                                       3716, 4868, 4913, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6278, 0, 3, 3788,
                                                                       3833, 4958, 5013, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6344, 0, 3, 3833,
                                                                       3878, 5013, 5068, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6410, 0, 3, 3878,
                                                                       3923, 5068, 5123, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6476, 0, 3, 3923,
                                                                       3968, 5123, 5178, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6542, 0, 3, 3968,
                                                                       4013, 5178, 5233, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 4013,
                                                                       4058, 5233, 5288, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6674, 0, 3, 4058,
                                                                       4103, 5288, 5343, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6740, 0, 3, 4103,
                                                                       4148, 5343, 5398, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6806, 0, 3, 4148,
                                                                       4193, 5398, 5453, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6872, 0, 3, 4193,
                                                                       4238, 5453, 5508, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6938, 0, 3, 4238,
                                                                       4283, 5508, 5563, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 7004, 0, 3, 4373,
                                                                       4418, 5618, 5673, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 7070, 0, 3, 4418,
                                                                       4463, 5673, 5728, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 7136, 0, 3, 4463,
                                                                       4508, 5728, 5783, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 7202, 0, 3, 4508,
                                                                       4553, 5783, 5838, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 7268, 0, 3, 4553,
                                                                       4598, 5838, 5893, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 7334, 0, 3, 4598,
                                                                       4643, 5893, 5948, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 4643,
                                                                       4688, 5948, 6003, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 7466, 0, 3, 4688,
                                                                       4733, 6003, 6058, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 7532, 0, 3, 4733,
                                                                       4778, 6058, 6113, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 7598, 0, 3, 4778,
                                                                       4823, 6113, 6168, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 7664, 0, 3, 4823,
                                                                       4868, 6168, 6223, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7730, 0, 3, 4958,
                                                                       5013, 6278, 6344, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7808, 0, 3, 5013,
                                                                       5068, 6344, 6410, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7886, 0, 3, 5068,
                                                                       5123, 6410, 6476, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7964, 0, 3, 5123,
                                                                       5178, 6476, 6542, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8042, 0, 3, 5178,
                                                                       5233, 6542, 6608, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8120, 0, 3, 5233,
                                                                       5288, 6608, 6674, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8198, 0, 3, 5288,
                                                                       5343, 6674, 6740, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8276, 0, 3, 5343,
                                                                       5398, 6740, 6806, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8354, 0, 3, 5398,
                                                                       5453, 6806, 6872, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8432, 0, 3, 5453,
                                                                       5508, 6872, 6938, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8510, 0, 3, 5618,
                                                                       5673, 7004, 7070, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8588, 0, 3, 5673,
                                                                       5728, 7070, 7136, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8666, 0, 3, 5728,
                                                                       5783, 7136, 7202, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8744, 0, 3, 5783,
                                                                       5838, 7202, 7268, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8822, 0, 3, 5838,
                                                                       5893, 7268, 7334, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8900, 0, 3, 5893,
                                                                       5948, 7334, 7400, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 8978, 0, 3, 5948,
                                                                       6003, 7400, 7466, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 9056, 0, 3, 6003,
                                                                       6058, 7466, 7532, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 9134, 0, 3, 6058,
                                                                       6113, 7532, 7598, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 9212, 0, 3, 6113,
                                                                       6168, 7598, 7664, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 9290, 0, 3, 6278,
                                                                       6344, 7730, 7808, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 9381, 0, 3, 6344,
                                                                       6410, 7808, 7886, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 9472, 0, 3, 6410,
                                                                       6476, 7886, 7964, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 9563, 0, 3, 6476,
                                                                       6542, 7964, 8042, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 9654, 0, 3, 6542,
                                                                       6608, 8042, 8120, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 9745, 0, 3, 6608,
                                                                       6674, 8120, 8198, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 9836, 0, 3, 6674,
                                                                       6740, 8198, 8276, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 9927, 0, 3, 6740,
                                                                       6806, 8276, 8354, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 10018, 0, 3, 6806,
                                                                       6872, 8354, 8432, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 10109, 0, 3, 7004,
                                                                       7070, 8510, 8588, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 10200, 0, 3, 7070,
                                                                       7136, 8588, 8666, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 10291, 0, 3, 7136,
                                                                       7202, 8666, 8744, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 10382, 0, 3, 7202,
                                                                       7268, 8744, 8822, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 10473, 0, 3, 7268,
                                                                       7334, 8822, 8900, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 10564, 0, 3, 7334,
                                                                       7400, 8900, 8978, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 10655, 0, 3, 7400,
                                                                       7466, 8978, 9056, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 10746, 0, 3, 7466,
                                                                       7532, 9056, 9134, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 10837, 0, 3, 7532,
                                                                       7598, 9134, 9212, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10928, 3, 9,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10931, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10934, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10937, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10940, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10943, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10946, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10949, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10952, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10955, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10958, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10961, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10964, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10967, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10970, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10973, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10976, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10979, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10982, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10985, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10988, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10991, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10994, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 10997, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11000, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11003, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11006, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11009, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11012, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11015, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11018, 3, 42,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11021, 3, 43,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11024, 3, 44,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11027, 3, 45,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11030, 3, 46,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11033, 3, 47,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11036, 3, 48,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 11039, 3, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11042, 3, 9, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11051, 3, 10, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11060, 3, 11, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11069, 3, 12, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11078, 3, 13, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11087, 3, 14, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11096, 3, 15, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11105, 3, 16, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11114, 3, 17, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11123, 3, 18, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11132, 3, 19, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11141, 3, 20, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11150, 3, 21, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11159, 3, 22, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11168, 3, 23, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11177, 3, 24, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11186, 3, 25, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11195, 3, 26, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11204, 3, 31, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11213, 3, 32, 119,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11222, 3, 33, 122,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11231, 3, 34, 125,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11240, 3, 35, 128,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11249, 3, 36, 131,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11258, 3, 37, 134,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11267, 3, 38, 137,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11276, 3, 39, 140,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11285, 3, 40, 143,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11294, 3, 41, 146,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11303, 3, 42, 149,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11312, 3, 43, 152,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11321, 3, 44, 155,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11330, 3, 45, 158,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11339, 3, 46, 161,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11348, 3, 47, 164,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 11357, 3, 48, 167,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11366, 3, 56, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11384, 3, 59, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11402, 3, 62, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11420, 3, 65, 200,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11438, 3, 68, 206,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11456, 3, 71, 212,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11474, 3, 74, 218,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11492, 3, 77, 224,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11510, 3, 80, 230,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11528, 3, 83, 236,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11546, 3, 86, 242,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11564, 3, 89, 248,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11582, 3, 92, 254,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11600, 3, 95, 260,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11618, 3, 98, 266,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11636, 3, 101,
                                                                       272, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11654, 3, 104,
                                                                       278, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11672, 3, 116,
                                                                       296, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11690, 3, 119,
                                                                       302, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11708, 3, 122,
                                                                       308, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11726, 3, 125,
                                                                       314, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11744, 3, 128,
                                                                       320, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11762, 3, 131,
                                                                       326, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11780, 3, 134,
                                                                       332, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11798, 3, 137,
                                                                       338, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11816, 3, 140,
                                                                       344, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11834, 3, 143,
                                                                       350, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11852, 3, 146,
                                                                       356, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11870, 3, 149,
                                                                       362, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11888, 3, 152,
                                                                       368, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11906, 3, 155,
                                                                       374, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11924, 3, 158,
                                                                       380, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11942, 3, 161,
                                                                       386, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 11960, 3, 164,
                                                                       392, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 11978, 3, 182,
                                                                       418, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12008, 3, 188,
                                                                       428, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12038, 3, 194,
                                                                       438, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12068, 3, 200,
                                                                       448, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12098, 3, 206,
                                                                       458, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12128, 3, 212,
                                                                       468, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12158, 3, 218,
                                                                       478, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12188, 3, 224,
                                                                       488, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12218, 3, 230,
                                                                       498, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12248, 3, 236,
                                                                       508, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12278, 3, 242,
                                                                       518, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12308, 3, 248,
                                                                       528, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12338, 3, 254,
                                                                       538, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12368, 3, 260,
                                                                       548, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12398, 3, 266,
                                                                       558, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12428, 3, 272,
                                                                       568, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12458, 3, 296,
                                                                       598, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12488, 3, 302,
                                                                       608, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12518, 3, 308,
                                                                       618, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12548, 3, 314,
                                                                       628, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12578, 3, 320,
                                                                       638, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12608, 3, 326,
                                                                       648, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12638, 3, 332,
                                                                       658, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12668, 3, 338,
                                                                       668, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12698, 3, 344,
                                                                       678, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12728, 3, 350,
                                                                       688, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12758, 3, 356,
                                                                       698, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12788, 3, 362,
                                                                       708, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12818, 3, 368,
                                                                       718, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12848, 3, 374,
                                                                       728, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12878, 3, 380,
                                                                       738, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 12908, 3, 386,
                                                                       748, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12938, 3, 418,
                                                                       788, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12983, 3, 428,
                                                                       803, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13028, 3, 438,
                                                                       818, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13073, 3, 448,
                                                                       833, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13118, 3, 458,
                                                                       848, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13163, 3, 468,
                                                                       863, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13208, 3, 478,
                                                                       878, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13253, 3, 488,
                                                                       893, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13298, 3, 498,
                                                                       908, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13343, 3, 508,
                                                                       923, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13388, 3, 518,
                                                                       938, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13433, 3, 528,
                                                                       953, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13478, 3, 538,
                                                                       968, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13523, 3, 548,
                                                                       983, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13568, 3, 558,
                                                                       998, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13613, 3, 598,
                                                                       1043, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13658, 3, 608,
                                                                       1058, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13703, 3, 618,
                                                                       1073, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13748, 3, 628,
                                                                       1088, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13793, 3, 638,
                                                                       1103, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13838, 3, 648,
                                                                       1118, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13883, 3, 658,
                                                                       1133, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13928, 3, 668,
                                                                       1148, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 13973, 3, 678,
                                                                       1163, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 14018, 3, 688,
                                                                       1178, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 14063, 3, 698,
                                                                       1193, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 14108, 3, 708,
                                                                       1208, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 14153, 3, 718,
                                                                       1223, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 14198, 3, 728,
                                                                       1238, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 14243, 3, 738,
                                                                       1253, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14288, 3, 788,
                                                                       1310, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14351, 3, 803,
                                                                       1331, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14414, 3, 818,
                                                                       1352, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14477, 3, 833,
                                                                       1373, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14540, 3, 848,
                                                                       1394, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14603, 3, 863,
                                                                       1415, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14666, 3, 878,
                                                                       1436, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14729, 3, 893,
                                                                       1457, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14792, 3, 908,
                                                                       1478, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14855, 3, 923,
                                                                       1499, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14918, 3, 938,
                                                                       1520, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14981, 3, 953,
                                                                       1541, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15044, 3, 968,
                                                                       1562, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15107, 3, 983,
                                                                       1583, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15170, 3, 1043,
                                                                       1646, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15233, 3, 1058,
                                                                       1667, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15296, 3, 1073,
                                                                       1688, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15359, 3, 1088,
                                                                       1709, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15422, 3, 1103,
                                                                       1730, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15485, 3, 1118,
                                                                       1751, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15548, 3, 1133,
                                                                       1772, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15611, 3, 1148,
                                                                       1793, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15674, 3, 1163,
                                                                       1814, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15737, 3, 1178,
                                                                       1835, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15800, 3, 1193,
                                                                       1856, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15863, 3, 1208,
                                                                       1877, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15926, 3, 1223,
                                                                       1898, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 15989, 3, 1238,
                                                                       1919, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16052, 3, 1310,
                                                                       1996, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16136, 3, 1331,
                                                                       2024, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16220, 3, 1352,
                                                                       2052, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16304, 3, 1373,
                                                                       2080, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16388, 3, 1394,
                                                                       2108, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16472, 3, 1415,
                                                                       2136, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16556, 3, 1436,
                                                                       2164, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16640, 3, 1457,
                                                                       2192, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16724, 3, 1478,
                                                                       2220, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16808, 3, 1499,
                                                                       2248, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16892, 3, 1520,
                                                                       2276, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16976, 3, 1541,
                                                                       2304, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17060, 3, 1562,
                                                                       2332, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17144, 3, 1646,
                                                                       2416, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17228, 3, 1667,
                                                                       2444, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17312, 3, 1688,
                                                                       2472, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17396, 3, 1709,
                                                                       2500, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17480, 3, 1730,
                                                                       2528, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17564, 3, 1751,
                                                                       2556, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17648, 3, 1772,
                                                                       2584, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17732, 3, 1793,
                                                                       2612, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17816, 3, 1814,
                                                                       2640, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17900, 3, 1835,
                                                                       2668, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 17984, 3, 1856,
                                                                       2696, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 18068, 3, 1877,
                                                                       2724, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 18152, 3, 1898,
                                                                       2752, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18236, 3, 1996,
                                                                       2852, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18344, 3, 2024,
                                                                       2888, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18452, 3, 2052,
                                                                       2924, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18560, 3, 2080,
                                                                       2960, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18668, 3, 2108,
                                                                       2996, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18776, 3, 2136,
                                                                       3032, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18884, 3, 2164,
                                                                       3068, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18992, 3, 2192,
                                                                       3104, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 19100, 3, 2220,
                                                                       3140, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 19208, 3, 2248,
                                                                       3176, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 19316, 3, 2276,
                                                                       3212, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 19424, 3, 2304,
                                                                       3248, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 19532, 3, 2416,
                                                                       3356, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 19640, 3, 2444,
                                                                       3392, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 19748, 3, 2472,
                                                                       3428, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 19856, 3, 2500,
                                                                       3464, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 19964, 3, 2528,
                                                                       3500, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 20072, 3, 2556,
                                                                       3536, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 20180, 3, 2584,
                                                                       3572, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 20288, 3, 2612,
                                                                       3608, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 20396, 3, 2640,
                                                                       3644, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 20504, 3, 2668,
                                                                       3680, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 20612, 3, 2696,
                                                                       3716, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 20720, 3, 2724,
                                                                       3752, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 20828, 3, 2852,
                                                                       3878, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 20963, 3, 2888,
                                                                       3923, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21098, 3, 2924,
                                                                       3968, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21233, 3, 2960,
                                                                       4013, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21368, 3, 2996,
                                                                       4058, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21503, 3, 3032,
                                                                       4103, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21638, 3, 3068,
                                                                       4148, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21773, 3, 3104,
                                                                       4193, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21908, 3, 3140,
                                                                       4238, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 22043, 3, 3176,
                                                                       4283, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 22178, 3, 3212,
                                                                       4328, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 22313, 3, 3356,
                                                                       4463, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 22448, 3, 3392,
                                                                       4508, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 22583, 3, 3428,
                                                                       4553, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 22718, 3, 3464,
                                                                       4598, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 22853, 3, 3500,
                                                                       4643, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 22988, 3, 3536,
                                                                       4688, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 23123, 3, 3572,
                                                                       4733, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 23258, 3, 3608,
                                                                       4778, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 23393, 3, 3644,
                                                                       4823, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 23528, 3, 3680,
                                                                       4868, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 23663, 3, 3716,
                                                                       4913, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 23798, 3, 3878,
                                                                       5068, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 23963, 3, 3923,
                                                                       5123, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24128, 3, 3968,
                                                                       5178, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24293, 3, 4013,
                                                                       5233, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24458, 3, 4058,
                                                                       5288, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24623, 3, 4103,
                                                                       5343, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24788, 3, 4148,
                                                                       5398, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24953, 3, 4193,
                                                                       5453, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 25118, 3, 4238,
                                                                       5508, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 25283, 3, 4283,
                                                                       5563, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 25448, 3, 4463,
                                                                       5728, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 25613, 3, 4508,
                                                                       5783, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 25778, 3, 4553,
                                                                       5838, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 25943, 3, 4598,
                                                                       5893, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 26108, 3, 4643,
                                                                       5948, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 26273, 3, 4688,
                                                                       6003, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 26438, 3, 4733,
                                                                       6058, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 26603, 3, 4778,
                                                                       6113, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 26768, 3, 4823,
                                                                       6168, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 26933, 3, 4868,
                                                                       6223, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 27098, 3, 5068,
                                                                       6410, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 27296, 3, 5123,
                                                                       6476, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 27494, 3, 5178,
                                                                       6542, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 27692, 3, 5233,
                                                                       6608, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 27890, 3, 5288,
                                                                       6674, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 28088, 3, 5343,
                                                                       6740, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 28286, 3, 5398,
                                                                       6806, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 28484, 3, 5453,
                                                                       6872, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 28682, 3, 5508,
                                                                       6938, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 28880, 3, 5728,
                                                                       7136, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 29078, 3, 5783,
                                                                       7202, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 29276, 3, 5838,
                                                                       7268, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 29474, 3, 5893,
                                                                       7334, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 29672, 3, 5948,
                                                                       7400, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 29870, 3, 6003,
                                                                       7466, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 30068, 3, 6058,
                                                                       7532, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 30266, 3, 6113,
                                                                       7598, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 30464, 3, 6168,
                                                                       7664, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 30662, 3, 6410,
                                                                       7886, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 30896, 3, 6476,
                                                                       7964, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 31130, 3, 6542,
                                                                       8042, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 31364, 3, 6608,
                                                                       8120, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 31598, 3, 6674,
                                                                       8198, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 31832, 3, 6740,
                                                                       8276, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 32066, 3, 6806,
                                                                       8354, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 32300, 3, 6872,
                                                                       8432, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 32534, 3, 7136,
                                                                       8666, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 32768, 3, 7202,
                                                                       8744, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 33002, 3, 7268,
                                                                       8822, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 33236, 3, 7334,
                                                                       8900, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 33470, 3, 7400,
                                                                       8978, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 33704, 3, 7466,
                                                                       9056, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 33938, 3, 7532,
                                                                       9134, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 34172, 3, 7598,
                                                                       9212, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 34406, 3, 7886,
                                                                       9472, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 34679, 3, 7964,
                                                                       9563, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 34952, 3, 8042,
                                                                       9654, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 35225, 3, 8120,
                                                                       9745, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 35498, 3, 8198,
                                                                       9836, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 35771, 3, 8276,
                                                                       9927, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 36044, 3, 8354,
                                                                       10018, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 36317, 3, 8666,
                                                                       10291, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 36590, 3, 8744,
                                                                       10382, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 36863, 3, 8822,
                                                                       10473, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 37136, 3, 8900,
                                                                       10564, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 37409, 3, 8978,
                                                                       10655, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 37682, 3, 9056,
                                                                       10746, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 37955, 3, 9134,
                                                                       10837, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38228, 3, 7, 8,
                                                                       10928, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38234, 3, 8, 9,
                                                                       10931, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38240, 3, 9, 10,
                                                                       10934, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38246, 3, 10, 11,
                                                                       10937, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38252, 3, 11, 12,
                                                                       10940, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38258, 3, 12, 13,
                                                                       10943, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38264, 3, 13, 14,
                                                                       10946, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38270, 3, 14, 15,
                                                                       10949, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38276, 3, 15, 16,
                                                                       10952, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38282, 3, 16, 17,
                                                                       10955, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38288, 3, 17, 18,
                                                                       10958, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38294, 3, 18, 19,
                                                                       10961, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38300, 3, 19, 20,
                                                                       10964, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38306, 3, 20, 21,
                                                                       10967, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38312, 3, 21, 22,
                                                                       10970, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38318, 3, 22, 23,
                                                                       10973, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38324, 3, 23, 24,
                                                                       10976, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38330, 3, 24, 25,
                                                                       10979, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38336, 3, 25, 26,
                                                                       10982, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38342, 3, 29, 30,
                                                                       10985, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38348, 3, 30, 31,
                                                                       10988, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38354, 3, 31, 32,
                                                                       10991, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38360, 3, 32, 33,
                                                                       10994, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38366, 3, 33, 34,
                                                                       10997, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38372, 3, 34, 35,
                                                                       11000, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38378, 3, 35, 36,
                                                                       11003, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38384, 3, 36, 37,
                                                                       11006, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38390, 3, 37, 38,
                                                                       11009, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38396, 3, 38, 39,
                                                                       11012, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38402, 3, 39, 40,
                                                                       11015, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38408, 3, 40, 41,
                                                                       11018, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38414, 3, 41, 42,
                                                                       11021, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38420, 3, 42, 43,
                                                                       11024, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38426, 3, 43, 44,
                                                                       11027, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38432, 3, 44, 45,
                                                                       11030, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38438, 3, 45, 46,
                                                                       11033, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38444, 3, 46, 47,
                                                                       11036, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 38450, 3, 47, 48,
                                                                       11039, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38456, 0, 3,
                                                                       38228, 10928, 38234,
                                                                       11042, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38474, 0, 3,
                                                                       38234, 10931, 38240,
                                                                       11051, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38492, 0, 3,
                                                                       38240, 10934, 38246,
                                                                       11060, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38510, 0, 3,
                                                                       38246, 10937, 38252,
                                                                       11069, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38528, 0, 3,
                                                                       38252, 10940, 38258,
                                                                       11078, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38546, 0, 3,
                                                                       38258, 10943, 38264,
                                                                       11087, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38564, 0, 3,
                                                                       38264, 10946, 38270,
                                                                       11096, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38582, 0, 3,
                                                                       38270, 10949, 38276,
                                                                       11105, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38600, 0, 3,
                                                                       38276, 10952, 38282,
                                                                       11114, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38618, 0, 3,
                                                                       38282, 10955, 38288,
                                                                       11123, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38636, 0, 3,
                                                                       38288, 10958, 38294,
                                                                       11132, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38654, 0, 3,
                                                                       38294, 10961, 38300,
                                                                       11141, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38672, 0, 3,
                                                                       38300, 10964, 38306,
                                                                       11150, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38690, 0, 3,
                                                                       38306, 10967, 38312,
                                                                       11159, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38708, 0, 3,
                                                                       38312, 10970, 38318,
                                                                       11168, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38726, 0, 3,
                                                                       38318, 10973, 38324,
                                                                       11177, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38744, 0, 3,
                                                                       38324, 10976, 38330,
                                                                       11186, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38762, 0, 3,
                                                                       38330, 10979, 38336,
                                                                       11195, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38780, 0, 3,
                                                                       38342, 10985, 38348,
                                                                       11204, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38798, 0, 3,
                                                                       38348, 10988, 38354,
                                                                       11213, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38816, 0, 3,
                                                                       38354, 10991, 38360,
                                                                       11222, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38834, 0, 3,
                                                                       38360, 10994, 38366,
                                                                       11231, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38852, 0, 3,
                                                                       38366, 10997, 38372,
                                                                       11240, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38870, 0, 3,
                                                                       38372, 11000, 38378,
                                                                       11249, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38888, 0, 3,
                                                                       38378, 11003, 38384,
                                                                       11258, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38906, 0, 3,
                                                                       38384, 11006, 38390,
                                                                       11267, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38924, 0, 3,
                                                                       38390, 11009, 38396,
                                                                       11276, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38942, 0, 3,
                                                                       38396, 11012, 38402,
                                                                       11285, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38960, 0, 3,
                                                                       38402, 11015, 38408,
                                                                       11294, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38978, 0, 3,
                                                                       38408, 11018, 38414,
                                                                       11303, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 38996, 0, 3,
                                                                       38414, 11021, 38420,
                                                                       11312, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 39014, 0, 3,
                                                                       38420, 11024, 38426,
                                                                       11321, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 39032, 0, 3,
                                                                       38426, 11027, 38432,
                                                                       11330, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 39050, 0, 3,
                                                                       38432, 11030, 38438,
                                                                       11339, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 39068, 0, 3,
                                                                       38438, 11033, 38444,
                                                                       11348, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 39086, 0, 3,
                                                                       38444, 11036, 38450,
                                                                       11357, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39104, 0, 3,
                                                                       38456, 11042, 38474, 170,
                                                                       176, 11366, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39140, 0, 3,
                                                                       38474, 11051, 38492, 176,
                                                                       182, 11384, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39176, 0, 3,
                                                                       38492, 11060, 38510, 182,
                                                                       188, 11402, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39212, 0, 3,
                                                                       38510, 11069, 38528, 188,
                                                                       194, 11420, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39248, 0, 3,
                                                                       38528, 11078, 38546, 194,
                                                                       200, 11438, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39284, 0, 3,
                                                                       38546, 11087, 38564, 200,
                                                                       206, 11456, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39320, 0, 3,
                                                                       38564, 11096, 38582, 206,
                                                                       212, 11474, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39356, 0, 3,
                                                                       38582, 11105, 38600, 212,
                                                                       218, 11492, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39392, 0, 3,
                                                                       38600, 11114, 38618, 218,
                                                                       224, 11510, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39428, 0, 3,
                                                                       38618, 11123, 38636, 224,
                                                                       230, 11528, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39464, 0, 3,
                                                                       38636, 11132, 38654, 230,
                                                                       236, 11546, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39500, 0, 3,
                                                                       38654, 11141, 38672, 236,
                                                                       242, 11564, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39536, 0, 3,
                                                                       38672, 11150, 38690, 242,
                                                                       248, 11582, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39572, 0, 3,
                                                                       38690, 11159, 38708, 248,
                                                                       254, 11600, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39608, 0, 3,
                                                                       38708, 11168, 38726, 254,
                                                                       260, 11618, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39644, 0, 3,
                                                                       38726, 11177, 38744, 260,
                                                                       266, 11636, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39680, 0, 3,
                                                                       38744, 11186, 38762, 266,
                                                                       272, 11654, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39716, 0, 3,
                                                                       38780, 11204, 38798, 284,
                                                                       290, 11672, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39752, 0, 3,
                                                                       38798, 11213, 38816, 290,
                                                                       296, 11690, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39788, 0, 3,
                                                                       38816, 11222, 38834, 296,
                                                                       302, 11708, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39824, 0, 3,
                                                                       38834, 11231, 38852, 302,
                                                                       308, 11726, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39860, 0, 3,
                                                                       38852, 11240, 38870, 308,
                                                                       314, 11744, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39896, 0, 3,
                                                                       38870, 11249, 38888, 314,
                                                                       320, 11762, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39932, 0, 3,
                                                                       38888, 11258, 38906, 320,
                                                                       326, 11780, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 39968, 0, 3,
                                                                       38906, 11267, 38924, 326,
                                                                       332, 11798, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 40004, 0, 3,
                                                                       38924, 11276, 38942, 332,
                                                                       338, 11816, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 40040, 0, 3,
                                                                       38942, 11285, 38960, 338,
                                                                       344, 11834, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 40076, 0, 3,
                                                                       38960, 11294, 38978, 344,
                                                                       350, 11852, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 40112, 0, 3,
                                                                       38978, 11303, 38996, 350,
                                                                       356, 11870, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 40148, 0, 3,
                                                                       38996, 11312, 39014, 356,
                                                                       362, 11888, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 40184, 0, 3,
                                                                       39014, 11321, 39032, 362,
                                                                       368, 11906, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 40220, 0, 3,
                                                                       39032, 11330, 39050, 368,
                                                                       374, 11924, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 40256, 0, 3,
                                                                       39050, 11339, 39068, 374,
                                                                       380, 11942, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 40292, 0, 3,
                                                                       39068, 11348, 39086, 380,
                                                                       386, 11960, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40328, 0, 3,
                                                                       39104, 11366, 39140, 398,
                                                                       408, 11978, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40388, 0, 3,
                                                                       39140, 11384, 39176, 408,
                                                                       418, 12008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40448, 0, 3,
                                                                       39176, 11402, 39212, 418,
                                                                       428, 12038, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40508, 0, 3,
                                                                       39212, 11420, 39248, 428,
                                                                       438, 12068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40568, 0, 3,
                                                                       39248, 11438, 39284, 438,
                                                                       448, 12098, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40628, 0, 3,
                                                                       39284, 11456, 39320, 448,
                                                                       458, 12128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40688, 0, 3,
                                                                       39320, 11474, 39356, 458,
                                                                       468, 12158, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40748, 0, 3,
                                                                       39356, 11492, 39392, 468,
                                                                       478, 12188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40808, 0, 3,
                                                                       39392, 11510, 39428, 478,
                                                                       488, 12218, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40868, 0, 3,
                                                                       39428, 11528, 39464, 488,
                                                                       498, 12248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40928, 0, 3,
                                                                       39464, 11546, 39500, 498,
                                                                       508, 12278, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 40988, 0, 3,
                                                                       39500, 11564, 39536, 508,
                                                                       518, 12308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41048, 0, 3,
                                                                       39536, 11582, 39572, 518,
                                                                       528, 12338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41108, 0, 3,
                                                                       39572, 11600, 39608, 528,
                                                                       538, 12368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41168, 0, 3,
                                                                       39608, 11618, 39644, 538,
                                                                       548, 12398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41228, 0, 3,
                                                                       39644, 11636, 39680, 548,
                                                                       558, 12428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41288, 0, 3,
                                                                       39716, 11672, 39752, 578,
                                                                       588, 12458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41348, 0, 3,
                                                                       39752, 11690, 39788, 588,
                                                                       598, 12488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41408, 0, 3,
                                                                       39788, 11708, 39824, 598,
                                                                       608, 12518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41468, 0, 3,
                                                                       39824, 11726, 39860, 608,
                                                                       618, 12548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41528, 0, 3,
                                                                       39860, 11744, 39896, 618,
                                                                       628, 12578, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41588, 0, 3,
                                                                       39896, 11762, 39932, 628,
                                                                       638, 12608, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41648, 0, 3,
                                                                       39932, 11780, 39968, 638,
                                                                       648, 12638, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41708, 0, 3,
                                                                       39968, 11798, 40004, 648,
                                                                       658, 12668, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41768, 0, 3,
                                                                       40004, 11816, 40040, 658,
                                                                       668, 12698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41828, 0, 3,
                                                                       40040, 11834, 40076, 668,
                                                                       678, 12728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41888, 0, 3,
                                                                       40076, 11852, 40112, 678,
                                                                       688, 12758, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 41948, 0, 3,
                                                                       40112, 11870, 40148, 688,
                                                                       698, 12788, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 42008, 0, 3,
                                                                       40148, 11888, 40184, 698,
                                                                       708, 12818, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 42068, 0, 3,
                                                                       40184, 11906, 40220, 708,
                                                                       718, 12848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 42128, 0, 3,
                                                                       40220, 11924, 40256, 718,
                                                                       728, 12878, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 42188, 0, 3,
                                                                       40256, 11942, 40292, 728,
                                                                       738, 12908, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42248, 0, 3,
                                                                       40328, 11978, 40388, 758,
                                                                       773, 12938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42338, 0, 3,
                                                                       40388, 12008, 40448, 773,
                                                                       788, 12983, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42428, 0, 3,
                                                                       40448, 12038, 40508, 788,
                                                                       803, 13028, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42518, 0, 3,
                                                                       40508, 12068, 40568, 803,
                                                                       818, 13073, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42608, 0, 3,
                                                                       40568, 12098, 40628, 818,
                                                                       833, 13118, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42698, 0, 3,
                                                                       40628, 12128, 40688, 833,
                                                                       848, 13163, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42788, 0, 3,
                                                                       40688, 12158, 40748, 848,
                                                                       863, 13208, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42878, 0, 3,
                                                                       40748, 12188, 40808, 863,
                                                                       878, 13253, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42968, 0, 3,
                                                                       40808, 12218, 40868, 878,
                                                                       893, 13298, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 43058, 0, 3,
                                                                       40868, 12248, 40928, 893,
                                                                       908, 13343, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 43148, 0, 3,
                                                                       40928, 12278, 40988, 908,
                                                                       923, 13388, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 43238, 0, 3,
                                                                       40988, 12308, 41048, 923,
                                                                       938, 13433, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 43328, 0, 3,
                                                                       41048, 12338, 41108, 938,
                                                                       953, 13478, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 43418, 0, 3,
                                                                       41108, 12368, 41168, 953,
                                                                       968, 13523, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 43508, 0, 3,
                                                                       41168, 12398, 41228, 968,
                                                                       983, 13568, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 43598, 0, 3,
                                                                       41288, 12458, 41348, 1013,
                                                                       1028, 13613, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 43688, 0, 3,
                                                                       41348, 12488, 41408, 1028,
                                                                       1043, 13658, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 43778, 0, 3,
                                                                       41408, 12518, 41468, 1043,
                                                                       1058, 13703, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 43868, 0, 3,
                                                                       41468, 12548, 41528, 1058,
                                                                       1073, 13748, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 43958, 0, 3,
                                                                       41528, 12578, 41588, 1073,
                                                                       1088, 13793, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 44048, 0, 3,
                                                                       41588, 12608, 41648, 1088,
                                                                       1103, 13838, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 44138, 0, 3,
                                                                       41648, 12638, 41708, 1103,
                                                                       1118, 13883, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 44228, 0, 3,
                                                                       41708, 12668, 41768, 1118,
                                                                       1133, 13928, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 44318, 0, 3,
                                                                       41768, 12698, 41828, 1133,
                                                                       1148, 13973, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 44408, 0, 3,
                                                                       41828, 12728, 41888, 1148,
                                                                       1163, 14018, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 44498, 0, 3,
                                                                       41888, 12758, 41948, 1163,
                                                                       1178, 14063, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 44588, 0, 3,
                                                                       41948, 12788, 42008, 1178,
                                                                       1193, 14108, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 44678, 0, 3,
                                                                       42008, 12818, 42068, 1193,
                                                                       1208, 14153, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 44768, 0, 3,
                                                                       42068, 12848, 42128, 1208,
                                                                       1223, 14198, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 44858, 0, 3,
                                                                       42128, 12878, 42188, 1223,
                                                                       1238, 14243, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 44948, 0, 3,
                                                                       42248, 12938, 42338, 1268,
                                                                       1289, 14288, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 45074, 0, 3,
                                                                       42338, 12983, 42428, 1289,
                                                                       1310, 14351, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 45200, 0, 3,
                                                                       42428, 13028, 42518, 1310,
                                                                       1331, 14414, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 45326, 0, 3,
                                                                       42518, 13073, 42608, 1331,
                                                                       1352, 14477, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 45452, 0, 3,
                                                                       42608, 13118, 42698, 1352,
                                                                       1373, 14540, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 45578, 0, 3,
                                                                       42698, 13163, 42788, 1373,
                                                                       1394, 14603, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 45704, 0, 3,
                                                                       42788, 13208, 42878, 1394,
                                                                       1415, 14666, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 45830, 0, 3,
                                                                       42878, 13253, 42968, 1415,
                                                                       1436, 14729, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 45956, 0, 3,
                                                                       42968, 13298, 43058, 1436,
                                                                       1457, 14792, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 46082, 0, 3,
                                                                       43058, 13343, 43148, 1457,
                                                                       1478, 14855, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 46208, 0, 3,
                                                                       43148, 13388, 43238, 1478,
                                                                       1499, 14918, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 46334, 0, 3,
                                                                       43238, 13433, 43328, 1499,
                                                                       1520, 14981, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 46460, 0, 3,
                                                                       43328, 13478, 43418, 1520,
                                                                       1541, 15044, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 46586, 0, 3,
                                                                       43418, 13523, 43508, 1541,
                                                                       1562, 15107, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 46712, 0, 3,
                                                                       43598, 13613, 43688, 1604,
                                                                       1625, 15170, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 46838, 0, 3,
                                                                       43688, 13658, 43778, 1625,
                                                                       1646, 15233, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 46964, 0, 3,
                                                                       43778, 13703, 43868, 1646,
                                                                       1667, 15296, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 47090, 0, 3,
                                                                       43868, 13748, 43958, 1667,
                                                                       1688, 15359, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 47216, 0, 3,
                                                                       43958, 13793, 44048, 1688,
                                                                       1709, 15422, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 47342, 0, 3,
                                                                       44048, 13838, 44138, 1709,
                                                                       1730, 15485, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 47468, 0, 3,
                                                                       44138, 13883, 44228, 1730,
                                                                       1751, 15548, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 47594, 0, 3,
                                                                       44228, 13928, 44318, 1751,
                                                                       1772, 15611, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 47720, 0, 3,
                                                                       44318, 13973, 44408, 1772,
                                                                       1793, 15674, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 47846, 0, 3,
                                                                       44408, 14018, 44498, 1793,
                                                                       1814, 15737, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 47972, 0, 3,
                                                                       44498, 14063, 44588, 1814,
                                                                       1835, 15800, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 48098, 0, 3,
                                                                       44588, 14108, 44678, 1835,
                                                                       1856, 15863, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 48224, 0, 3,
                                                                       44678, 14153, 44768, 1856,
                                                                       1877, 15926, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 48350, 0, 3,
                                                                       44768, 14198, 44858, 1877,
                                                                       1898, 15989, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 48476, 0, 3,
                                                                       44948, 14288, 45074, 1940,
                                                                       1968, 16052, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 48644, 0, 3,
                                                                       45074, 14351, 45200, 1968,
                                                                       1996, 16136, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 48812, 0, 3,
                                                                       45200, 14414, 45326, 1996,
                                                                       2024, 16220, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 48980, 0, 3,
                                                                       45326, 14477, 45452, 2024,
                                                                       2052, 16304, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 49148, 0, 3,
                                                                       45452, 14540, 45578, 2052,
                                                                       2080, 16388, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 49316, 0, 3,
                                                                       45578, 14603, 45704, 2080,
                                                                       2108, 16472, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 49484, 0, 3,
                                                                       45704, 14666, 45830, 2108,
                                                                       2136, 16556, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 49652, 0, 3,
                                                                       45830, 14729, 45956, 2136,
                                                                       2164, 16640, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 49820, 0, 3,
                                                                       45956, 14792, 46082, 2164,
                                                                       2192, 16724, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 49988, 0, 3,
                                                                       46082, 14855, 46208, 2192,
                                                                       2220, 16808, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 50156, 0, 3,
                                                                       46208, 14918, 46334, 2220,
                                                                       2248, 16892, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 50324, 0, 3,
                                                                       46334, 14981, 46460, 2248,
                                                                       2276, 16976, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 50492, 0, 3,
                                                                       46460, 15044, 46586, 2276,
                                                                       2304, 17060, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 50660, 0, 3,
                                                                       46712, 15170, 46838, 2360,
                                                                       2388, 17144, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 50828, 0, 3,
                                                                       46838, 15233, 46964, 2388,
                                                                       2416, 17228, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 50996, 0, 3,
                                                                       46964, 15296, 47090, 2416,
                                                                       2444, 17312, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 51164, 0, 3,
                                                                       47090, 15359, 47216, 2444,
                                                                       2472, 17396, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 51332, 0, 3,
                                                                       47216, 15422, 47342, 2472,
                                                                       2500, 17480, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 51500, 0, 3,
                                                                       47342, 15485, 47468, 2500,
                                                                       2528, 17564, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 51668, 0, 3,
                                                                       47468, 15548, 47594, 2528,
                                                                       2556, 17648, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 51836, 0, 3,
                                                                       47594, 15611, 47720, 2556,
                                                                       2584, 17732, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 52004, 0, 3,
                                                                       47720, 15674, 47846, 2584,
                                                                       2612, 17816, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 52172, 0, 3,
                                                                       47846, 15737, 47972, 2612,
                                                                       2640, 17900, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 52340, 0, 3,
                                                                       47972, 15800, 48098, 2640,
                                                                       2668, 17984, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 52508, 0, 3,
                                                                       48098, 15863, 48224, 2668,
                                                                       2696, 18068, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 52676, 0, 3,
                                                                       48224, 15926, 48350, 2696,
                                                                       2724, 18152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 52844, 0, 3,
                                                                       48476, 16052, 48644, 2780,
                                                                       2816, 18236, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 53060, 0, 3,
                                                                       48644, 16136, 48812, 2816,
                                                                       2852, 18344, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 53276, 0, 3,
                                                                       48812, 16220, 48980, 2852,
                                                                       2888, 18452, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 53492, 0, 3,
                                                                       48980, 16304, 49148, 2888,
                                                                       2924, 18560, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 53708, 0, 3,
                                                                       49148, 16388, 49316, 2924,
                                                                       2960, 18668, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 53924, 0, 3,
                                                                       49316, 16472, 49484, 2960,
                                                                       2996, 18776, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 54140, 0, 3,
                                                                       49484, 16556, 49652, 2996,
                                                                       3032, 18884, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 54356, 0, 3,
                                                                       49652, 16640, 49820, 3032,
                                                                       3068, 18992, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 54572, 0, 3,
                                                                       49820, 16724, 49988, 3068,
                                                                       3104, 19100, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 54788, 0, 3,
                                                                       49988, 16808, 50156, 3104,
                                                                       3140, 19208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 55004, 0, 3,
                                                                       50156, 16892, 50324, 3140,
                                                                       3176, 19316, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 55220, 0, 3,
                                                                       50324, 16976, 50492, 3176,
                                                                       3212, 19424, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 55436, 0, 3,
                                                                       50660, 17144, 50828, 3284,
                                                                       3320, 19532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 55652, 0, 3,
                                                                       50828, 17228, 50996, 3320,
                                                                       3356, 19640, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 55868, 0, 3,
                                                                       50996, 17312, 51164, 3356,
                                                                       3392, 19748, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 56084, 0, 3,
                                                                       51164, 17396, 51332, 3392,
                                                                       3428, 19856, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 56300, 0, 3,
                                                                       51332, 17480, 51500, 3428,
                                                                       3464, 19964, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 56516, 0, 3,
                                                                       51500, 17564, 51668, 3464,
                                                                       3500, 20072, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 56732, 0, 3,
                                                                       51668, 17648, 51836, 3500,
                                                                       3536, 20180, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 56948, 0, 3,
                                                                       51836, 17732, 52004, 3536,
                                                                       3572, 20288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 57164, 0, 3,
                                                                       52004, 17816, 52172, 3572,
                                                                       3608, 20396, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 57380, 0, 3,
                                                                       52172, 17900, 52340, 3608,
                                                                       3644, 20504, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 57596, 0, 3,
                                                                       52340, 17984, 52508, 3644,
                                                                       3680, 20612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 57812, 0, 3,
                                                                       52508, 18068, 52676, 3680,
                                                                       3716, 20720, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 58028, 0, 3,
                                                                       52844, 18236, 53060, 3788,
                                                                       3833, 20828, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 58298, 0, 3,
                                                                       53060, 18344, 53276, 3833,
                                                                       3878, 20963, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 58568, 0, 3,
                                                                       53276, 18452, 53492, 3878,
                                                                       3923, 21098, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 58838, 0, 3,
                                                                       53492, 18560, 53708, 3923,
                                                                       3968, 21233, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 59108, 0, 3,
                                                                       53708, 18668, 53924, 3968,
                                                                       4013, 21368, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 59378, 0, 3,
                                                                       53924, 18776, 54140, 4013,
                                                                       4058, 21503, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 59648, 0, 3,
                                                                       54140, 18884, 54356, 4058,
                                                                       4103, 21638, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 59918, 0, 3,
                                                                       54356, 18992, 54572, 4103,
                                                                       4148, 21773, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 60188, 0, 3,
                                                                       54572, 19100, 54788, 4148,
                                                                       4193, 21908, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 60458, 0, 3,
                                                                       54788, 19208, 55004, 4193,
                                                                       4238, 22043, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 60728, 0, 3,
                                                                       55004, 19316, 55220, 4238,
                                                                       4283, 22178, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 60998, 0, 3,
                                                                       55436, 19532, 55652, 4373,
                                                                       4418, 22313, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 61268, 0, 3,
                                                                       55652, 19640, 55868, 4418,
                                                                       4463, 22448, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 61538, 0, 3,
                                                                       55868, 19748, 56084, 4463,
                                                                       4508, 22583, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 61808, 0, 3,
                                                                       56084, 19856, 56300, 4508,
                                                                       4553, 22718, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 62078, 0, 3,
                                                                       56300, 19964, 56516, 4553,
                                                                       4598, 22853, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 62348, 0, 3,
                                                                       56516, 20072, 56732, 4598,
                                                                       4643, 22988, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 62618, 0, 3,
                                                                       56732, 20180, 56948, 4643,
                                                                       4688, 23123, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 62888, 0, 3,
                                                                       56948, 20288, 57164, 4688,
                                                                       4733, 23258, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 63158, 0, 3,
                                                                       57164, 20396, 57380, 4733,
                                                                       4778, 23393, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 63428, 0, 3,
                                                                       57380, 20504, 57596, 4778,
                                                                       4823, 23528, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 63698, 0, 3,
                                                                       57596, 20612, 57812, 4823,
                                                                       4868, 23663, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 63968, 0, 3,
                                                                       58028, 20828, 58298, 4958,
                                                                       5013, 23798, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 64298, 0, 3,
                                                                       58298, 20963, 58568, 5013,
                                                                       5068, 23963, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 64628, 0, 3,
                                                                       58568, 21098, 58838, 5068,
                                                                       5123, 24128, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 64958, 0, 3,
                                                                       58838, 21233, 59108, 5123,
                                                                       5178, 24293, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 65288, 0, 3,
                                                                       59108, 21368, 59378, 5178,
                                                                       5233, 24458, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 65618, 0, 3,
                                                                       59378, 21503, 59648, 5233,
                                                                       5288, 24623, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 65948, 0, 3,
                                                                       59648, 21638, 59918, 5288,
                                                                       5343, 24788, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 66278, 0, 3,
                                                                       59918, 21773, 60188, 5343,
                                                                       5398, 24953, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 66608, 0, 3,
                                                                       60188, 21908, 60458, 5398,
                                                                       5453, 25118, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 66938, 0, 3,
                                                                       60458, 22043, 60728, 5453,
                                                                       5508, 25283, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 67268, 0, 3,
                                                                       60998, 22313, 61268, 5618,
                                                                       5673, 25448, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 67598, 0, 3,
                                                                       61268, 22448, 61538, 5673,
                                                                       5728, 25613, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 67928, 0, 3,
                                                                       61538, 22583, 61808, 5728,
                                                                       5783, 25778, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 68258, 0, 3,
                                                                       61808, 22718, 62078, 5783,
                                                                       5838, 25943, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 68588, 0, 3,
                                                                       62078, 22853, 62348, 5838,
                                                                       5893, 26108, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 68918, 0, 3,
                                                                       62348, 22988, 62618, 5893,
                                                                       5948, 26273, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 69248, 0, 3,
                                                                       62618, 23123, 62888, 5948,
                                                                       6003, 26438, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 69578, 0, 3,
                                                                       62888, 23258, 63158, 6003,
                                                                       6058, 26603, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 69908, 0, 3,
                                                                       63158, 23393, 63428, 6058,
                                                                       6113, 26768, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 70238, 0, 3,
                                                                       63428, 23528, 63698, 6113,
                                                                       6168, 26933, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 70568, 0, 3,
                                                                       63968, 23798, 64298, 6278,
                                                                       6344, 27098, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 70964, 0, 3,
                                                                       64298, 23963, 64628, 6344,
                                                                       6410, 27296, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 71360, 0, 3,
                                                                       64628, 24128, 64958, 6410,
                                                                       6476, 27494, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 71756, 0, 3,
                                                                       64958, 24293, 65288, 6476,
                                                                       6542, 27692, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 72152, 0, 3,
                                                                       65288, 24458, 65618, 6542,
                                                                       6608, 27890, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 72548, 0, 3,
                                                                       65618, 24623, 65948, 6608,
                                                                       6674, 28088, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 72944, 0, 3,
                                                                       65948, 24788, 66278, 6674,
                                                                       6740, 28286, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 73340, 0, 3,
                                                                       66278, 24953, 66608, 6740,
                                                                       6806, 28484, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 73736, 0, 3,
                                                                       66608, 25118, 66938, 6806,
                                                                       6872, 28682, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 74132, 0, 3,
                                                                       67268, 25448, 67598, 7004,
                                                                       7070, 28880, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 74528, 0, 3,
                                                                       67598, 25613, 67928, 7070,
                                                                       7136, 29078, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 74924, 0, 3,
                                                                       67928, 25778, 68258, 7136,
                                                                       7202, 29276, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 75320, 0, 3,
                                                                       68258, 25943, 68588, 7202,
                                                                       7268, 29474, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 75716, 0, 3,
                                                                       68588, 26108, 68918, 7268,
                                                                       7334, 29672, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 76112, 0, 3,
                                                                       68918, 26273, 69248, 7334,
                                                                       7400, 29870, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 76508, 0, 3,
                                                                       69248, 26438, 69578, 7400,
                                                                       7466, 30068, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 76904, 0, 3,
                                                                       69578, 26603, 69908, 7466,
                                                                       7532, 30266, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 77300, 0, 3,
                                                                       69908, 26768, 70238, 7532,
                                                                       7598, 30464, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 77696, 0, 3,
                                                                       70568, 27098, 70964, 7730,
                                                                       7808, 30662, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 78164, 0, 3,
                                                                       70964, 27296, 71360, 7808,
                                                                       7886, 30896, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 78632, 0, 3,
                                                                       71360, 27494, 71756, 7886,
                                                                       7964, 31130, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 79100, 0, 3,
                                                                       71756, 27692, 72152, 7964,
                                                                       8042, 31364, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 79568, 0, 3,
                                                                       72152, 27890, 72548, 8042,
                                                                       8120, 31598, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 80036, 0, 3,
                                                                       72548, 28088, 72944, 8120,
                                                                       8198, 31832, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 80504, 0, 3,
                                                                       72944, 28286, 73340, 8198,
                                                                       8276, 32066, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 80972, 0, 3,
                                                                       73340, 28484, 73736, 8276,
                                                                       8354, 32300, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 81440, 0, 3,
                                                                       74132, 28880, 74528, 8510,
                                                                       8588, 32534, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 81908, 0, 3,
                                                                       74528, 29078, 74924, 8588,
                                                                       8666, 32768, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 82376, 0, 3,
                                                                       74924, 29276, 75320, 8666,
                                                                       8744, 33002, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 82844, 0, 3,
                                                                       75320, 29474, 75716, 8744,
                                                                       8822, 33236, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 83312, 0, 3,
                                                                       75716, 29672, 76112, 8822,
                                                                       8900, 33470, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 83780, 0, 3,
                                                                       76112, 29870, 76508, 8900,
                                                                       8978, 33704, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 84248, 0, 3,
                                                                       76508, 30068, 76904, 8978,
                                                                       9056, 33938, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 84716, 0, 3,
                                                                       76904, 30266, 77300, 9056,
                                                                       9134, 34172, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 85184, 0, 3,
                                                                       77696, 30662, 78164, 9290,
                                                                       9381, 34406, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 85730, 0, 3,
                                                                       78164, 30896, 78632, 9381,
                                                                       9472, 34679, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 86276, 0, 3,
                                                                       78632, 31130, 79100, 9472,
                                                                       9563, 34952, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 86822, 0, 3,
                                                                       79100, 31364, 79568, 9563,
                                                                       9654, 35225, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 87368, 0, 3,
                                                                       79568, 31598, 80036, 9654,
                                                                       9745, 35498, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 87914, 0, 3,
                                                                       80036, 31832, 80504, 9745,
                                                                       9836, 35771, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 88460, 0, 3,
                                                                       80504, 32066, 80972, 9836,
                                                                       9927, 36044, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 89006, 0, 3,
                                                                       81440, 32534, 81908,
                                                                       10109, 10200, 36317,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 89552, 0, 3,
                                                                       81908, 32768, 82376,
                                                                       10200, 10291, 36590,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 90098, 0, 3,
                                                                       82376, 33002, 82844,
                                                                       10291, 10382, 36863,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 90644, 0, 3,
                                                                       82844, 33236, 83312,
                                                                       10382, 10473, 37136,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 91190, 0, 3,
                                                                       83312, 33470, 83780,
                                                                       10473, 10564, 37409,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 91736, 0, 3,
                                                                       83780, 33704, 84248,
                                                                       10564, 10655, 37682,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 92282, 0, 3,
                                                                       84248, 33938, 84716,
                                                                       10655, 10746, 37955,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92828, 3, 10928,
                                                                       10931, 38240, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92838, 3, 10931,
                                                                       10934, 38246, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92848, 3, 10934,
                                                                       10937, 38252, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92858, 3, 10937,
                                                                       10940, 38258, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92868, 3, 10940,
                                                                       10943, 38264, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92878, 3, 10943,
                                                                       10946, 38270, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92888, 3, 10946,
                                                                       10949, 38276, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92898, 3, 10949,
                                                                       10952, 38282, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92908, 3, 10952,
                                                                       10955, 38288, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92918, 3, 10955,
                                                                       10958, 38294, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92928, 3, 10958,
                                                                       10961, 38300, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92938, 3, 10961,
                                                                       10964, 38306, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92948, 3, 10964,
                                                                       10967, 38312, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92958, 3, 10967,
                                                                       10970, 38318, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92968, 3, 10970,
                                                                       10973, 38324, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92978, 3, 10973,
                                                                       10976, 38330, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92988, 3, 10976,
                                                                       10979, 38336, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 92998, 3, 10985,
                                                                       10988, 38354, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93008, 3, 10988,
                                                                       10991, 38360, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93018, 3, 10991,
                                                                       10994, 38366, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93028, 3, 10994,
                                                                       10997, 38372, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93038, 3, 10997,
                                                                       11000, 38378, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93048, 3, 11000,
                                                                       11003, 38384, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93058, 3, 11003,
                                                                       11006, 38390, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93068, 3, 11006,
                                                                       11009, 38396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93078, 3, 11009,
                                                                       11012, 38402, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93088, 3, 11012,
                                                                       11015, 38408, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93098, 3, 11015,
                                                                       11018, 38414, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93108, 3, 11018,
                                                                       11021, 38420, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93118, 3, 11021,
                                                                       11024, 38426, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93128, 3, 11024,
                                                                       11027, 38432, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93138, 3, 11027,
                                                                       11030, 38438, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93148, 3, 11030,
                                                                       11033, 38444, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 93158, 3, 11033,
                                                                       11036, 38450, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93168, 0, 3,
                                                                       92828, 38240, 92838,
                                                                       11042, 11051, 38492,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93198, 0, 3,
                                                                       92838, 38246, 92848,
                                                                       11051, 11060, 38510,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93228, 0, 3,
                                                                       92848, 38252, 92858,
                                                                       11060, 11069, 38528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93258, 0, 3,
                                                                       92858, 38258, 92868,
                                                                       11069, 11078, 38546,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93288, 0, 3,
                                                                       92868, 38264, 92878,
                                                                       11078, 11087, 38564,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93318, 0, 3,
                                                                       92878, 38270, 92888,
                                                                       11087, 11096, 38582,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93348, 0, 3,
                                                                       92888, 38276, 92898,
                                                                       11096, 11105, 38600,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93378, 0, 3,
                                                                       92898, 38282, 92908,
                                                                       11105, 11114, 38618,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93408, 0, 3,
                                                                       92908, 38288, 92918,
                                                                       11114, 11123, 38636,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93438, 0, 3,
                                                                       92918, 38294, 92928,
                                                                       11123, 11132, 38654,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93468, 0, 3,
                                                                       92928, 38300, 92938,
                                                                       11132, 11141, 38672,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93498, 0, 3,
                                                                       92938, 38306, 92948,
                                                                       11141, 11150, 38690,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93528, 0, 3,
                                                                       92948, 38312, 92958,
                                                                       11150, 11159, 38708,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93558, 0, 3,
                                                                       92958, 38318, 92968,
                                                                       11159, 11168, 38726,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93588, 0, 3,
                                                                       92968, 38324, 92978,
                                                                       11168, 11177, 38744,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93618, 0, 3,
                                                                       92978, 38330, 92988,
                                                                       11177, 11186, 38762,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93648, 0, 3,
                                                                       92998, 38354, 93008,
                                                                       11204, 11213, 38816,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93678, 0, 3,
                                                                       93008, 38360, 93018,
                                                                       11213, 11222, 38834,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93708, 0, 3,
                                                                       93018, 38366, 93028,
                                                                       11222, 11231, 38852,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93738, 0, 3,
                                                                       93028, 38372, 93038,
                                                                       11231, 11240, 38870,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93768, 0, 3,
                                                                       93038, 38378, 93048,
                                                                       11240, 11249, 38888,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93798, 0, 3,
                                                                       93048, 38384, 93058,
                                                                       11249, 11258, 38906,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93828, 0, 3,
                                                                       93058, 38390, 93068,
                                                                       11258, 11267, 38924,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93858, 0, 3,
                                                                       93068, 38396, 93078,
                                                                       11267, 11276, 38942,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93888, 0, 3,
                                                                       93078, 38402, 93088,
                                                                       11276, 11285, 38960,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93918, 0, 3,
                                                                       93088, 38408, 93098,
                                                                       11285, 11294, 38978,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93948, 0, 3,
                                                                       93098, 38414, 93108,
                                                                       11294, 11303, 38996,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 93978, 0, 3,
                                                                       93108, 38420, 93118,
                                                                       11303, 11312, 39014,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 94008, 0, 3,
                                                                       93118, 38426, 93128,
                                                                       11312, 11321, 39032,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 94038, 0, 3,
                                                                       93128, 38432, 93138,
                                                                       11321, 11330, 39050,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 94068, 0, 3,
                                                                       93138, 38438, 93148,
                                                                       11330, 11339, 39068,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 94098, 0, 3,
                                                                       93148, 38444, 93158,
                                                                       11339, 11348, 39086,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94128, 0, 3,
                                                                       93168, 38492, 93198,
                                                                       11366, 11384, 39176,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94188, 0, 3,
                                                                       93198, 38510, 93228,
                                                                       11384, 11402, 39212,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94248, 0, 3,
                                                                       93228, 38528, 93258,
                                                                       11402, 11420, 39248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94308, 0, 3,
                                                                       93258, 38546, 93288,
                                                                       11420, 11438, 39284,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94368, 0, 3,
                                                                       93288, 38564, 93318,
                                                                       11438, 11456, 39320,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94428, 0, 3,
                                                                       93318, 38582, 93348,
                                                                       11456, 11474, 39356,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94488, 0, 3,
                                                                       93348, 38600, 93378,
                                                                       11474, 11492, 39392,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94548, 0, 3,
                                                                       93378, 38618, 93408,
                                                                       11492, 11510, 39428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94608, 0, 3,
                                                                       93408, 38636, 93438,
                                                                       11510, 11528, 39464,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94668, 0, 3,
                                                                       93438, 38654, 93468,
                                                                       11528, 11546, 39500,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94728, 0, 3,
                                                                       93468, 38672, 93498,
                                                                       11546, 11564, 39536,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94788, 0, 3,
                                                                       93498, 38690, 93528,
                                                                       11564, 11582, 39572,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94848, 0, 3,
                                                                       93528, 38708, 93558,
                                                                       11582, 11600, 39608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94908, 0, 3,
                                                                       93558, 38726, 93588,
                                                                       11600, 11618, 39644,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 94968, 0, 3,
                                                                       93588, 38744, 93618,
                                                                       11618, 11636, 39680,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95028, 0, 3,
                                                                       93648, 38816, 93678,
                                                                       11672, 11690, 39788,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95088, 0, 3,
                                                                       93678, 38834, 93708,
                                                                       11690, 11708, 39824,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95148, 0, 3,
                                                                       93708, 38852, 93738,
                                                                       11708, 11726, 39860,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95208, 0, 3,
                                                                       93738, 38870, 93768,
                                                                       11726, 11744, 39896,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95268, 0, 3,
                                                                       93768, 38888, 93798,
                                                                       11744, 11762, 39932,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95328, 0, 3,
                                                                       93798, 38906, 93828,
                                                                       11762, 11780, 39968,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95388, 0, 3,
                                                                       93828, 38924, 93858,
                                                                       11780, 11798, 40004,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95448, 0, 3,
                                                                       93858, 38942, 93888,
                                                                       11798, 11816, 40040,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95508, 0, 3,
                                                                       93888, 38960, 93918,
                                                                       11816, 11834, 40076,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95568, 0, 3,
                                                                       93918, 38978, 93948,
                                                                       11834, 11852, 40112,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95628, 0, 3,
                                                                       93948, 38996, 93978,
                                                                       11852, 11870, 40148,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95688, 0, 3,
                                                                       93978, 39014, 94008,
                                                                       11870, 11888, 40184,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95748, 0, 3,
                                                                       94008, 39032, 94038,
                                                                       11888, 11906, 40220,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95808, 0, 3,
                                                                       94038, 39050, 94068,
                                                                       11906, 11924, 40256,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 95868, 0, 3,
                                                                       94068, 39068, 94098,
                                                                       11924, 11942, 40292,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 95928, 0, 3,
                                                                       94128, 39176, 94188,
                                                                       11978, 12008, 40448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 96028, 0, 3,
                                                                       94188, 39212, 94248,
                                                                       12008, 12038, 40508,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 96128, 0, 3,
                                                                       94248, 39248, 94308,
                                                                       12038, 12068, 40568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 96228, 0, 3,
                                                                       94308, 39284, 94368,
                                                                       12068, 12098, 40628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 96328, 0, 3,
                                                                       94368, 39320, 94428,
                                                                       12098, 12128, 40688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 96428, 0, 3,
                                                                       94428, 39356, 94488,
                                                                       12128, 12158, 40748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 96528, 0, 3,
                                                                       94488, 39392, 94548,
                                                                       12158, 12188, 40808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 96628, 0, 3,
                                                                       94548, 39428, 94608,
                                                                       12188, 12218, 40868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 96728, 0, 3,
                                                                       94608, 39464, 94668,
                                                                       12218, 12248, 40928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 96828, 0, 3,
                                                                       94668, 39500, 94728,
                                                                       12248, 12278, 40988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 96928, 0, 3,
                                                                       94728, 39536, 94788,
                                                                       12278, 12308, 41048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 97028, 0, 3,
                                                                       94788, 39572, 94848,
                                                                       12308, 12338, 41108,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 97128, 0, 3,
                                                                       94848, 39608, 94908,
                                                                       12338, 12368, 41168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 97228, 0, 3,
                                                                       94908, 39644, 94968,
                                                                       12368, 12398, 41228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 97328, 0, 3,
                                                                       95028, 39788, 95088,
                                                                       12458, 12488, 41408,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 97428, 0, 3,
                                                                       95088, 39824, 95148,
                                                                       12488, 12518, 41468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 97528, 0, 3,
                                                                       95148, 39860, 95208,
                                                                       12518, 12548, 41528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 97628, 0, 3,
                                                                       95208, 39896, 95268,
                                                                       12548, 12578, 41588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 97728, 0, 3,
                                                                       95268, 39932, 95328,
                                                                       12578, 12608, 41648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 97828, 0, 3,
                                                                       95328, 39968, 95388,
                                                                       12608, 12638, 41708,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 97928, 0, 3,
                                                                       95388, 40004, 95448,
                                                                       12638, 12668, 41768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 98028, 0, 3,
                                                                       95448, 40040, 95508,
                                                                       12668, 12698, 41828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 98128, 0, 3,
                                                                       95508, 40076, 95568,
                                                                       12698, 12728, 41888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 98228, 0, 3,
                                                                       95568, 40112, 95628,
                                                                       12728, 12758, 41948,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 98328, 0, 3,
                                                                       95628, 40148, 95688,
                                                                       12758, 12788, 42008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 98428, 0, 3,
                                                                       95688, 40184, 95748,
                                                                       12788, 12818, 42068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 98528, 0, 3,
                                                                       95748, 40220, 95808,
                                                                       12818, 12848, 42128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 98628, 0, 3,
                                                                       95808, 40256, 95868,
                                                                       12848, 12878, 42188,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 98728, 0, 3,
                                                                       95928, 40448, 96028,
                                                                       12938, 12983, 42428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 98878, 0, 3,
                                                                       96028, 40508, 96128,
                                                                       12983, 13028, 42518,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 99028, 0, 3,
                                                                       96128, 40568, 96228,
                                                                       13028, 13073, 42608,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 99178, 0, 3,
                                                                       96228, 40628, 96328,
                                                                       13073, 13118, 42698,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 99328, 0, 3,
                                                                       96328, 40688, 96428,
                                                                       13118, 13163, 42788,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 99478, 0, 3,
                                                                       96428, 40748, 96528,
                                                                       13163, 13208, 42878,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 99628, 0, 3,
                                                                       96528, 40808, 96628,
                                                                       13208, 13253, 42968,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 99778, 0, 3,
                                                                       96628, 40868, 96728,
                                                                       13253, 13298, 43058,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 99928, 0, 3,
                                                                       96728, 40928, 96828,
                                                                       13298, 13343, 43148,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 100078, 0, 3,
                                                                       96828, 40988, 96928,
                                                                       13343, 13388, 43238,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 100228, 0, 3,
                                                                       96928, 41048, 97028,
                                                                       13388, 13433, 43328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 100378, 0, 3,
                                                                       97028, 41108, 97128,
                                                                       13433, 13478, 43418,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 100528, 0, 3,
                                                                       97128, 41168, 97228,
                                                                       13478, 13523, 43508,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 100678, 0, 3,
                                                                       97328, 41408, 97428,
                                                                       13613, 13658, 43778,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 100828, 0, 3,
                                                                       97428, 41468, 97528,
                                                                       13658, 13703, 43868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 100978, 0, 3,
                                                                       97528, 41528, 97628,
                                                                       13703, 13748, 43958,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 101128, 0, 3,
                                                                       97628, 41588, 97728,
                                                                       13748, 13793, 44048,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 101278, 0, 3,
                                                                       97728, 41648, 97828,
                                                                       13793, 13838, 44138,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 101428, 0, 3,
                                                                       97828, 41708, 97928,
                                                                       13838, 13883, 44228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 101578, 0, 3,
                                                                       97928, 41768, 98028,
                                                                       13883, 13928, 44318,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 101728, 0, 3,
                                                                       98028, 41828, 98128,
                                                                       13928, 13973, 44408,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 101878, 0, 3,
                                                                       98128, 41888, 98228,
                                                                       13973, 14018, 44498,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 102028, 0, 3,
                                                                       98228, 41948, 98328,
                                                                       14018, 14063, 44588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 102178, 0, 3,
                                                                       98328, 42008, 98428,
                                                                       14063, 14108, 44678,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 102328, 0, 3,
                                                                       98428, 42068, 98528,
                                                                       14108, 14153, 44768,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 102478, 0, 3,
                                                                       98528, 42128, 98628,
                                                                       14153, 14198, 44858,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 102628, 0, 3,
                                                                       98728, 42428, 98878,
                                                                       14288, 14351, 45200,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 102838, 0, 3,
                                                                       98878, 42518, 99028,
                                                                       14351, 14414, 45326,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 103048, 0, 3,
                                                                       99028, 42608, 99178,
                                                                       14414, 14477, 45452,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 103258, 0, 3,
                                                                       99178, 42698, 99328,
                                                                       14477, 14540, 45578,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 103468, 0, 3,
                                                                       99328, 42788, 99478,
                                                                       14540, 14603, 45704,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 103678, 0, 3,
                                                                       99478, 42878, 99628,
                                                                       14603, 14666, 45830,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 103888, 0, 3,
                                                                       99628, 42968, 99778,
                                                                       14666, 14729, 45956,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 104098, 0, 3,
                                                                       99778, 43058, 99928,
                                                                       14729, 14792, 46082,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 104308, 0, 3,
                                                                       99928, 43148, 100078,
                                                                       14792, 14855, 46208,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 104518, 0, 3,
                                                                       100078, 43238, 100228,
                                                                       14855, 14918, 46334,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 104728, 0, 3,
                                                                       100228, 43328, 100378,
                                                                       14918, 14981, 46460,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 104938, 0, 3,
                                                                       100378, 43418, 100528,
                                                                       14981, 15044, 46586,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 105148, 0, 3,
                                                                       100678, 43778, 100828,
                                                                       15170, 15233, 46964,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 105358, 0, 3,
                                                                       100828, 43868, 100978,
                                                                       15233, 15296, 47090,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 105568, 0, 3,
                                                                       100978, 43958, 101128,
                                                                       15296, 15359, 47216,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 105778, 0, 3,
                                                                       101128, 44048, 101278,
                                                                       15359, 15422, 47342,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 105988, 0, 3,
                                                                       101278, 44138, 101428,
                                                                       15422, 15485, 47468,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 106198, 0, 3,
                                                                       101428, 44228, 101578,
                                                                       15485, 15548, 47594,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 106408, 0, 3,
                                                                       101578, 44318, 101728,
                                                                       15548, 15611, 47720,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 106618, 0, 3,
                                                                       101728, 44408, 101878,
                                                                       15611, 15674, 47846,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 106828, 0, 3,
                                                                       101878, 44498, 102028,
                                                                       15674, 15737, 47972,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 107038, 0, 3,
                                                                       102028, 44588, 102178,
                                                                       15737, 15800, 48098,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 107248, 0, 3,
                                                                       102178, 44678, 102328,
                                                                       15800, 15863, 48224,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 107458, 0, 3,
                                                                       102328, 44768, 102478,
                                                                       15863, 15926, 48350,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 107668, 0, 3,
                                                                       102628, 45200, 102838,
                                                                       16052, 16136, 48812,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 107948, 0, 3,
                                                                       102838, 45326, 103048,
                                                                       16136, 16220, 48980,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 108228, 0, 3,
                                                                       103048, 45452, 103258,
                                                                       16220, 16304, 49148,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 108508, 0, 3,
                                                                       103258, 45578, 103468,
                                                                       16304, 16388, 49316,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 108788, 0, 3,
                                                                       103468, 45704, 103678,
                                                                       16388, 16472, 49484,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 109068, 0, 3,
                                                                       103678, 45830, 103888,
                                                                       16472, 16556, 49652,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 109348, 0, 3,
                                                                       103888, 45956, 104098,
                                                                       16556, 16640, 49820,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 109628, 0, 3,
                                                                       104098, 46082, 104308,
                                                                       16640, 16724, 49988,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 109908, 0, 3,
                                                                       104308, 46208, 104518,
                                                                       16724, 16808, 50156,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 110188, 0, 3,
                                                                       104518, 46334, 104728,
                                                                       16808, 16892, 50324,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 110468, 0, 3,
                                                                       104728, 46460, 104938,
                                                                       16892, 16976, 50492,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 110748, 0, 3,
                                                                       105148, 46964, 105358,
                                                                       17144, 17228, 50996,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 111028, 0, 3,
                                                                       105358, 47090, 105568,
                                                                       17228, 17312, 51164,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 111308, 0, 3,
                                                                       105568, 47216, 105778,
                                                                       17312, 17396, 51332,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 111588, 0, 3,
                                                                       105778, 47342, 105988,
                                                                       17396, 17480, 51500,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 111868, 0, 3,
                                                                       105988, 47468, 106198,
                                                                       17480, 17564, 51668,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 112148, 0, 3,
                                                                       106198, 47594, 106408,
                                                                       17564, 17648, 51836,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 112428, 0, 3,
                                                                       106408, 47720, 106618,
                                                                       17648, 17732, 52004,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 112708, 0, 3,
                                                                       106618, 47846, 106828,
                                                                       17732, 17816, 52172,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 112988, 0, 3,
                                                                       106828, 47972, 107038,
                                                                       17816, 17900, 52340,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 113268, 0, 3,
                                                                       107038, 48098, 107248,
                                                                       17900, 17984, 52508,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 113548, 0, 3,
                                                                       107248, 48224, 107458,
                                                                       17984, 18068, 52676,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 113828, 0, 3,
                                                                       107668, 48812, 107948,
                                                                       18236, 18344, 53276,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 114188, 0, 3,
                                                                       107948, 48980, 108228,
                                                                       18344, 18452, 53492,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 114548, 0, 3,
                                                                       108228, 49148, 108508,
                                                                       18452, 18560, 53708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 114908, 0, 3,
                                                                       108508, 49316, 108788,
                                                                       18560, 18668, 53924,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 115268, 0, 3,
                                                                       108788, 49484, 109068,
                                                                       18668, 18776, 54140,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 115628, 0, 3,
                                                                       109068, 49652, 109348,
                                                                       18776, 18884, 54356,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 115988, 0, 3,
                                                                       109348, 49820, 109628,
                                                                       18884, 18992, 54572,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 116348, 0, 3,
                                                                       109628, 49988, 109908,
                                                                       18992, 19100, 54788,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 116708, 0, 3,
                                                                       109908, 50156, 110188,
                                                                       19100, 19208, 55004,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 117068, 0, 3,
                                                                       110188, 50324, 110468,
                                                                       19208, 19316, 55220,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 117428, 0, 3,
                                                                       110748, 50996, 111028,
                                                                       19532, 19640, 55868,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 117788, 0, 3,
                                                                       111028, 51164, 111308,
                                                                       19640, 19748, 56084,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 118148, 0, 3,
                                                                       111308, 51332, 111588,
                                                                       19748, 19856, 56300,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 118508, 0, 3,
                                                                       111588, 51500, 111868,
                                                                       19856, 19964, 56516,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 118868, 0, 3,
                                                                       111868, 51668, 112148,
                                                                       19964, 20072, 56732,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 119228, 0, 3,
                                                                       112148, 51836, 112428,
                                                                       20072, 20180, 56948,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 119588, 0, 3,
                                                                       112428, 52004, 112708,
                                                                       20180, 20288, 57164,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 119948, 0, 3,
                                                                       112708, 52172, 112988,
                                                                       20288, 20396, 57380,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 120308, 0, 3,
                                                                       112988, 52340, 113268,
                                                                       20396, 20504, 57596,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 120668, 0, 3,
                                                                       113268, 52508, 113548,
                                                                       20504, 20612, 57812,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 121028, 0, 3,
                                                                       113828, 53276, 114188,
                                                                       20828, 20963, 58568,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 121478, 0, 3,
                                                                       114188, 53492, 114548,
                                                                       20963, 21098, 58838,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 121928, 0, 3,
                                                                       114548, 53708, 114908,
                                                                       21098, 21233, 59108,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 122378, 0, 3,
                                                                       114908, 53924, 115268,
                                                                       21233, 21368, 59378,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 122828, 0, 3,
                                                                       115268, 54140, 115628,
                                                                       21368, 21503, 59648,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 123278, 0, 3,
                                                                       115628, 54356, 115988,
                                                                       21503, 21638, 59918,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 123728, 0, 3,
                                                                       115988, 54572, 116348,
                                                                       21638, 21773, 60188,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 124178, 0, 3,
                                                                       116348, 54788, 116708,
                                                                       21773, 21908, 60458,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 124628, 0, 3,
                                                                       116708, 55004, 117068,
                                                                       21908, 22043, 60728,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 125078, 0, 3,
                                                                       117428, 55868, 117788,
                                                                       22313, 22448, 61538,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 125528, 0, 3,
                                                                       117788, 56084, 118148,
                                                                       22448, 22583, 61808,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 125978, 0, 3,
                                                                       118148, 56300, 118508,
                                                                       22583, 22718, 62078,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 126428, 0, 3,
                                                                       118508, 56516, 118868,
                                                                       22718, 22853, 62348,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 126878, 0, 3,
                                                                       118868, 56732, 119228,
                                                                       22853, 22988, 62618,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 127328, 0, 3,
                                                                       119228, 56948, 119588,
                                                                       22988, 23123, 62888,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 127778, 0, 3,
                                                                       119588, 57164, 119948,
                                                                       23123, 23258, 63158,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 128228, 0, 3,
                                                                       119948, 57380, 120308,
                                                                       23258, 23393, 63428,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 128678, 0, 3,
                                                                       120308, 57596, 120668,
                                                                       23393, 23528, 63698,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 129128, 0, 3,
                                                                       121028, 58568, 121478,
                                                                       23798, 23963, 64628,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 129678, 0, 3,
                                                                       121478, 58838, 121928,
                                                                       23963, 24128, 64958,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 130228, 0, 3,
                                                                       121928, 59108, 122378,
                                                                       24128, 24293, 65288,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 130778, 0, 3,
                                                                       122378, 59378, 122828,
                                                                       24293, 24458, 65618,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 131328, 0, 3,
                                                                       122828, 59648, 123278,
                                                                       24458, 24623, 65948,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 131878, 0, 3,
                                                                       123278, 59918, 123728,
                                                                       24623, 24788, 66278,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 132428, 0, 3,
                                                                       123728, 60188, 124178,
                                                                       24788, 24953, 66608,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 132978, 0, 3,
                                                                       124178, 60458, 124628,
                                                                       24953, 25118, 66938,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 133528, 0, 3,
                                                                       125078, 61538, 125528,
                                                                       25448, 25613, 67928,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 134078, 0, 3,
                                                                       125528, 61808, 125978,
                                                                       25613, 25778, 68258,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 134628, 0, 3,
                                                                       125978, 62078, 126428,
                                                                       25778, 25943, 68588,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 135178, 0, 3,
                                                                       126428, 62348, 126878,
                                                                       25943, 26108, 68918,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 135728, 0, 3,
                                                                       126878, 62618, 127328,
                                                                       26108, 26273, 69248,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 136278, 0, 3,
                                                                       127328, 62888, 127778,
                                                                       26273, 26438, 69578,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 136828, 0, 3,
                                                                       127778, 63158, 128228,
                                                                       26438, 26603, 69908,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 137378, 0, 3,
                                                                       128228, 63428, 128678,
                                                                       26603, 26768, 70238,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 137928, 0, 3,
                                                                       129128, 64628, 129678,
                                                                       27098, 27296, 71360,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 138588, 0, 3,
                                                                       129678, 64958, 130228,
                                                                       27296, 27494, 71756,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 139248, 0, 3,
                                                                       130228, 65288, 130778,
                                                                       27494, 27692, 72152,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 139908, 0, 3,
                                                                       130778, 65618, 131328,
                                                                       27692, 27890, 72548,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 140568, 0, 3,
                                                                       131328, 65948, 131878,
                                                                       27890, 28088, 72944,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 141228, 0, 3,
                                                                       131878, 66278, 132428,
                                                                       28088, 28286, 73340,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 141888, 0, 3,
                                                                       132428, 66608, 132978,
                                                                       28286, 28484, 73736,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 142548, 0, 3,
                                                                       133528, 67928, 134078,
                                                                       28880, 29078, 74924,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 143208, 0, 3,
                                                                       134078, 68258, 134628,
                                                                       29078, 29276, 75320,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 143868, 0, 3,
                                                                       134628, 68588, 135178,
                                                                       29276, 29474, 75716,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 144528, 0, 3,
                                                                       135178, 68918, 135728,
                                                                       29474, 29672, 76112,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 145188, 0, 3,
                                                                       135728, 69248, 136278,
                                                                       29672, 29870, 76508,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 145848, 0, 3,
                                                                       136278, 69578, 136828,
                                                                       29870, 30068, 76904,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 146508, 0, 3,
                                                                       136828, 69908, 137378,
                                                                       30068, 30266, 77300,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 147168, 0, 3,
                                                                       137928, 71360, 138588,
                                                                       30662, 30896, 78632,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 147948, 0, 3,
                                                                       138588, 71756, 139248,
                                                                       30896, 31130, 79100,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 148728, 0, 3,
                                                                       139248, 72152, 139908,
                                                                       31130, 31364, 79568,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 149508, 0, 3,
                                                                       139908, 72548, 140568,
                                                                       31364, 31598, 80036,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 150288, 0, 3,
                                                                       140568, 72944, 141228,
                                                                       31598, 31832, 80504,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 151068, 0, 3,
                                                                       141228, 73340, 141888,
                                                                       31832, 32066, 80972,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 151848, 0, 3,
                                                                       142548, 74924, 143208,
                                                                       32534, 32768, 82376,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 152628, 0, 3,
                                                                       143208, 75320, 143868,
                                                                       32768, 33002, 82844,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 153408, 0, 3,
                                                                       143868, 75716, 144528,
                                                                       33002, 33236, 83312,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 154188, 0, 3,
                                                                       144528, 76112, 145188,
                                                                       33236, 33470, 83780,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 154968, 0, 3,
                                                                       145188, 76508, 145848,
                                                                       33470, 33704, 84248,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 155748, 0, 3,
                                                                       145848, 76904, 146508,
                                                                       33704, 33938, 84716,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 156528, 0, 3,
                                                                       147168, 78632, 147948,
                                                                       34406, 34679, 86276,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 157438, 0, 3,
                                                                       147948, 79100, 148728,
                                                                       34679, 34952, 86822,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 158348, 0, 3,
                                                                       148728, 79568, 149508,
                                                                       34952, 35225, 87368,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 159258, 0, 3,
                                                                       149508, 80036, 150288,
                                                                       35225, 35498, 87914,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 160168, 0, 3,
                                                                       150288, 80504, 151068,
                                                                       35498, 35771, 88460,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 161078, 0, 3,
                                                                       151848, 82376, 152628,
                                                                       36317, 36590, 90098,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 161988, 0, 3,
                                                                       152628, 82844, 153408,
                                                                       36590, 36863, 90644,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 162898, 0, 3,
                                                                       153408, 83312, 154188,
                                                                       36863, 37136, 91190,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 163808, 0, 3,
                                                                       154188, 83780, 154968,
                                                                       37136, 37409, 91736,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 164718, 0, 3,
                                                                       154968, 84248, 155748,
                                                                       37409, 37682, 92282,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165628, 3, 38228,
                                                                       38234, 92828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165643, 3, 38234,
                                                                       38240, 92838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165658, 3, 38240,
                                                                       38246, 92848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165673, 3, 38246,
                                                                       38252, 92858, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165688, 3, 38252,
                                                                       38258, 92868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165703, 3, 38258,
                                                                       38264, 92878, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165718, 3, 38264,
                                                                       38270, 92888, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165733, 3, 38270,
                                                                       38276, 92898, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165748, 3, 38276,
                                                                       38282, 92908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165763, 3, 38282,
                                                                       38288, 92918, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165778, 3, 38288,
                                                                       38294, 92928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165793, 3, 38294,
                                                                       38300, 92938, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165808, 3, 38300,
                                                                       38306, 92948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165823, 3, 38306,
                                                                       38312, 92958, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165838, 3, 38312,
                                                                       38318, 92968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165853, 3, 38318,
                                                                       38324, 92978, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165868, 3, 38324,
                                                                       38330, 92988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165883, 3, 38342,
                                                                       38348, 92998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165898, 3, 38348,
                                                                       38354, 93008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165913, 3, 38354,
                                                                       38360, 93018, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165928, 3, 38360,
                                                                       38366, 93028, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165943, 3, 38366,
                                                                       38372, 93038, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165958, 3, 38372,
                                                                       38378, 93048, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165973, 3, 38378,
                                                                       38384, 93058, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 165988, 3, 38384,
                                                                       38390, 93068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 166003, 3, 38390,
                                                                       38396, 93078, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 166018, 3, 38396,
                                                                       38402, 93088, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 166033, 3, 38402,
                                                                       38408, 93098, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 166048, 3, 38408,
                                                                       38414, 93108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 166063, 3, 38414,
                                                                       38420, 93118, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 166078, 3, 38420,
                                                                       38426, 93128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 166093, 3, 38426,
                                                                       38432, 93138, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 166108, 3, 38432,
                                                                       38438, 93148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 166123, 3, 38438,
                                                                       38444, 93158, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166138, 0, 3,
                                                                       165628, 92828, 165643,
                                                                       38456, 38474, 93168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166183, 0, 3,
                                                                       165643, 92838, 165658,
                                                                       38474, 38492, 93198,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166228, 0, 3,
                                                                       165658, 92848, 165673,
                                                                       38492, 38510, 93228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166273, 0, 3,
                                                                       165673, 92858, 165688,
                                                                       38510, 38528, 93258,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166318, 0, 3,
                                                                       165688, 92868, 165703,
                                                                       38528, 38546, 93288,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166363, 0, 3,
                                                                       165703, 92878, 165718,
                                                                       38546, 38564, 93318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166408, 0, 3,
                                                                       165718, 92888, 165733,
                                                                       38564, 38582, 93348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166453, 0, 3,
                                                                       165733, 92898, 165748,
                                                                       38582, 38600, 93378,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166498, 0, 3,
                                                                       165748, 92908, 165763,
                                                                       38600, 38618, 93408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166543, 0, 3,
                                                                       165763, 92918, 165778,
                                                                       38618, 38636, 93438,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166588, 0, 3,
                                                                       165778, 92928, 165793,
                                                                       38636, 38654, 93468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166633, 0, 3,
                                                                       165793, 92938, 165808,
                                                                       38654, 38672, 93498,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166678, 0, 3,
                                                                       165808, 92948, 165823,
                                                                       38672, 38690, 93528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166723, 0, 3,
                                                                       165823, 92958, 165838,
                                                                       38690, 38708, 93558,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166768, 0, 3,
                                                                       165838, 92968, 165853,
                                                                       38708, 38726, 93588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166813, 0, 3,
                                                                       165853, 92978, 165868,
                                                                       38726, 38744, 93618,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166858, 0, 3,
                                                                       165883, 92998, 165898,
                                                                       38780, 38798, 93648,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166903, 0, 3,
                                                                       165898, 93008, 165913,
                                                                       38798, 38816, 93678,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166948, 0, 3,
                                                                       165913, 93018, 165928,
                                                                       38816, 38834, 93708,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 166993, 0, 3,
                                                                       165928, 93028, 165943,
                                                                       38834, 38852, 93738,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167038, 0, 3,
                                                                       165943, 93038, 165958,
                                                                       38852, 38870, 93768,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167083, 0, 3,
                                                                       165958, 93048, 165973,
                                                                       38870, 38888, 93798,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167128, 0, 3,
                                                                       165973, 93058, 165988,
                                                                       38888, 38906, 93828,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167173, 0, 3,
                                                                       165988, 93068, 166003,
                                                                       38906, 38924, 93858,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167218, 0, 3,
                                                                       166003, 93078, 166018,
                                                                       38924, 38942, 93888,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167263, 0, 3,
                                                                       166018, 93088, 166033,
                                                                       38942, 38960, 93918,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167308, 0, 3,
                                                                       166033, 93098, 166048,
                                                                       38960, 38978, 93948,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167353, 0, 3,
                                                                       166048, 93108, 166063,
                                                                       38978, 38996, 93978,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167398, 0, 3,
                                                                       166063, 93118, 166078,
                                                                       38996, 39014, 94008,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167443, 0, 3,
                                                                       166078, 93128, 166093,
                                                                       39014, 39032, 94038,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167488, 0, 3,
                                                                       166093, 93138, 166108,
                                                                       39032, 39050, 94068,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 167533, 0, 3,
                                                                       166108, 93148, 166123,
                                                                       39050, 39068, 94098,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 167578, 0, 3,
                                                                       166138, 93168, 166183,
                                                                       39104, 39140, 94128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 167668, 0, 3,
                                                                       166183, 93198, 166228,
                                                                       39140, 39176, 94188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 167758, 0, 3,
                                                                       166228, 93228, 166273,
                                                                       39176, 39212, 94248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 167848, 0, 3,
                                                                       166273, 93258, 166318,
                                                                       39212, 39248, 94308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 167938, 0, 3,
                                                                       166318, 93288, 166363,
                                                                       39248, 39284, 94368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 168028, 0, 3,
                                                                       166363, 93318, 166408,
                                                                       39284, 39320, 94428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 168118, 0, 3,
                                                                       166408, 93348, 166453,
                                                                       39320, 39356, 94488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 168208, 0, 3,
                                                                       166453, 93378, 166498,
                                                                       39356, 39392, 94548,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 168298, 0, 3,
                                                                       166498, 93408, 166543,
                                                                       39392, 39428, 94608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 168388, 0, 3,
                                                                       166543, 93438, 166588,
                                                                       39428, 39464, 94668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 168478, 0, 3,
                                                                       166588, 93468, 166633,
                                                                       39464, 39500, 94728,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 168568, 0, 3,
                                                                       166633, 93498, 166678,
                                                                       39500, 39536, 94788,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 168658, 0, 3,
                                                                       166678, 93528, 166723,
                                                                       39536, 39572, 94848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 168748, 0, 3,
                                                                       166723, 93558, 166768,
                                                                       39572, 39608, 94908,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 168838, 0, 3,
                                                                       166768, 93588, 166813,
                                                                       39608, 39644, 94968,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 168928, 0, 3,
                                                                       166858, 93648, 166903,
                                                                       39716, 39752, 95028,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 169018, 0, 3,
                                                                       166903, 93678, 166948,
                                                                       39752, 39788, 95088,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 169108, 0, 3,
                                                                       166948, 93708, 166993,
                                                                       39788, 39824, 95148,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 169198, 0, 3,
                                                                       166993, 93738, 167038,
                                                                       39824, 39860, 95208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 169288, 0, 3,
                                                                       167038, 93768, 167083,
                                                                       39860, 39896, 95268,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 169378, 0, 3,
                                                                       167083, 93798, 167128,
                                                                       39896, 39932, 95328,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 169468, 0, 3,
                                                                       167128, 93828, 167173,
                                                                       39932, 39968, 95388,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 169558, 0, 3,
                                                                       167173, 93858, 167218,
                                                                       39968, 40004, 95448,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 169648, 0, 3,
                                                                       167218, 93888, 167263,
                                                                       40004, 40040, 95508,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 169738, 0, 3,
                                                                       167263, 93918, 167308,
                                                                       40040, 40076, 95568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 169828, 0, 3,
                                                                       167308, 93948, 167353,
                                                                       40076, 40112, 95628,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 169918, 0, 3,
                                                                       167353, 93978, 167398,
                                                                       40112, 40148, 95688,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 170008, 0, 3,
                                                                       167398, 94008, 167443,
                                                                       40148, 40184, 95748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 170098, 0, 3,
                                                                       167443, 94038, 167488,
                                                                       40184, 40220, 95808,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 170188, 0, 3,
                                                                       167488, 94068, 167533,
                                                                       40220, 40256, 95868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 170278, 0, 3,
                                                                       167578, 94128, 167668,
                                                                       40328, 40388, 95928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 170428, 0, 3,
                                                                       167668, 94188, 167758,
                                                                       40388, 40448, 96028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 170578, 0, 3,
                                                                       167758, 94248, 167848,
                                                                       40448, 40508, 96128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 170728, 0, 3,
                                                                       167848, 94308, 167938,
                                                                       40508, 40568, 96228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 170878, 0, 3,
                                                                       167938, 94368, 168028,
                                                                       40568, 40628, 96328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 171028, 0, 3,
                                                                       168028, 94428, 168118,
                                                                       40628, 40688, 96428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 171178, 0, 3,
                                                                       168118, 94488, 168208,
                                                                       40688, 40748, 96528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 171328, 0, 3,
                                                                       168208, 94548, 168298,
                                                                       40748, 40808, 96628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 171478, 0, 3,
                                                                       168298, 94608, 168388,
                                                                       40808, 40868, 96728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 171628, 0, 3,
                                                                       168388, 94668, 168478,
                                                                       40868, 40928, 96828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 171778, 0, 3,
                                                                       168478, 94728, 168568,
                                                                       40928, 40988, 96928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 171928, 0, 3,
                                                                       168568, 94788, 168658,
                                                                       40988, 41048, 97028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 172078, 0, 3,
                                                                       168658, 94848, 168748,
                                                                       41048, 41108, 97128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 172228, 0, 3,
                                                                       168748, 94908, 168838,
                                                                       41108, 41168, 97228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 172378, 0, 3,
                                                                       168928, 95028, 169018,
                                                                       41288, 41348, 97328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 172528, 0, 3,
                                                                       169018, 95088, 169108,
                                                                       41348, 41408, 97428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 172678, 0, 3,
                                                                       169108, 95148, 169198,
                                                                       41408, 41468, 97528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 172828, 0, 3,
                                                                       169198, 95208, 169288,
                                                                       41468, 41528, 97628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 172978, 0, 3,
                                                                       169288, 95268, 169378,
                                                                       41528, 41588, 97728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 173128, 0, 3,
                                                                       169378, 95328, 169468,
                                                                       41588, 41648, 97828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 173278, 0, 3,
                                                                       169468, 95388, 169558,
                                                                       41648, 41708, 97928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 173428, 0, 3,
                                                                       169558, 95448, 169648,
                                                                       41708, 41768, 98028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 173578, 0, 3,
                                                                       169648, 95508, 169738,
                                                                       41768, 41828, 98128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 173728, 0, 3,
                                                                       169738, 95568, 169828,
                                                                       41828, 41888, 98228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 173878, 0, 3,
                                                                       169828, 95628, 169918,
                                                                       41888, 41948, 98328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 174028, 0, 3,
                                                                       169918, 95688, 170008,
                                                                       41948, 42008, 98428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 174178, 0, 3,
                                                                       170008, 95748, 170098,
                                                                       42008, 42068, 98528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 174328, 0, 3,
                                                                       170098, 95808, 170188,
                                                                       42068, 42128, 98628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 174478, 0, 3,
                                                                       170278, 95928, 170428,
                                                                       42248, 42338, 98728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 174703, 0, 3,
                                                                       170428, 96028, 170578,
                                                                       42338, 42428, 98878,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 174928, 0, 3,
                                                                       170578, 96128, 170728,
                                                                       42428, 42518, 99028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 175153, 0, 3,
                                                                       170728, 96228, 170878,
                                                                       42518, 42608, 99178,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 175378, 0, 3,
                                                                       170878, 96328, 171028,
                                                                       42608, 42698, 99328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 175603, 0, 3,
                                                                       171028, 96428, 171178,
                                                                       42698, 42788, 99478,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 175828, 0, 3,
                                                                       171178, 96528, 171328,
                                                                       42788, 42878, 99628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 176053, 0, 3,
                                                                       171328, 96628, 171478,
                                                                       42878, 42968, 99778,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 176278, 0, 3,
                                                                       171478, 96728, 171628,
                                                                       42968, 43058, 99928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 176503, 0, 3,
                                                                       171628, 96828, 171778,
                                                                       43058, 43148, 100078,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 176728, 0, 3,
                                                                       171778, 96928, 171928,
                                                                       43148, 43238, 100228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 176953, 0, 3,
                                                                       171928, 97028, 172078,
                                                                       43238, 43328, 100378,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 177178, 0, 3,
                                                                       172078, 97128, 172228,
                                                                       43328, 43418, 100528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 177403, 0, 3,
                                                                       172378, 97328, 172528,
                                                                       43598, 43688, 100678,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 177628, 0, 3,
                                                                       172528, 97428, 172678,
                                                                       43688, 43778, 100828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 177853, 0, 3,
                                                                       172678, 97528, 172828,
                                                                       43778, 43868, 100978,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 178078, 0, 3,
                                                                       172828, 97628, 172978,
                                                                       43868, 43958, 101128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 178303, 0, 3,
                                                                       172978, 97728, 173128,
                                                                       43958, 44048, 101278,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 178528, 0, 3,
                                                                       173128, 97828, 173278,
                                                                       44048, 44138, 101428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 178753, 0, 3,
                                                                       173278, 97928, 173428,
                                                                       44138, 44228, 101578,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 178978, 0, 3,
                                                                       173428, 98028, 173578,
                                                                       44228, 44318, 101728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 179203, 0, 3,
                                                                       173578, 98128, 173728,
                                                                       44318, 44408, 101878,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 179428, 0, 3,
                                                                       173728, 98228, 173878,
                                                                       44408, 44498, 102028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 179653, 0, 3,
                                                                       173878, 98328, 174028,
                                                                       44498, 44588, 102178,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 179878, 0, 3,
                                                                       174028, 98428, 174178,
                                                                       44588, 44678, 102328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 180103, 0, 3,
                                                                       174178, 98528, 174328,
                                                                       44678, 44768, 102478,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 180328, 0, 3,
                                                                       174478, 98728, 174703,
                                                                       44948, 45074, 102628,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 180643, 0, 3,
                                                                       174703, 98878, 174928,
                                                                       45074, 45200, 102838,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 180958, 0, 3,
                                                                       174928, 99028, 175153,
                                                                       45200, 45326, 103048,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 181273, 0, 3,
                                                                       175153, 99178, 175378,
                                                                       45326, 45452, 103258,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 181588, 0, 3,
                                                                       175378, 99328, 175603,
                                                                       45452, 45578, 103468,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 181903, 0, 3,
                                                                       175603, 99478, 175828,
                                                                       45578, 45704, 103678,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 182218, 0, 3,
                                                                       175828, 99628, 176053,
                                                                       45704, 45830, 103888,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 182533, 0, 3,
                                                                       176053, 99778, 176278,
                                                                       45830, 45956, 104098,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 182848, 0, 3,
                                                                       176278, 99928, 176503,
                                                                       45956, 46082, 104308,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 183163, 0, 3,
                                                                       176503, 100078, 176728,
                                                                       46082, 46208, 104518,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 183478, 0, 3,
                                                                       176728, 100228, 176953,
                                                                       46208, 46334, 104728,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 183793, 0, 3,
                                                                       176953, 100378, 177178,
                                                                       46334, 46460, 104938,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 184108, 0, 3,
                                                                       177403, 100678, 177628,
                                                                       46712, 46838, 105148,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 184423, 0, 3,
                                                                       177628, 100828, 177853,
                                                                       46838, 46964, 105358,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 184738, 0, 3,
                                                                       177853, 100978, 178078,
                                                                       46964, 47090, 105568,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 185053, 0, 3,
                                                                       178078, 101128, 178303,
                                                                       47090, 47216, 105778,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 185368, 0, 3,
                                                                       178303, 101278, 178528,
                                                                       47216, 47342, 105988,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 185683, 0, 3,
                                                                       178528, 101428, 178753,
                                                                       47342, 47468, 106198,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 185998, 0, 3,
                                                                       178753, 101578, 178978,
                                                                       47468, 47594, 106408,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 186313, 0, 3,
                                                                       178978, 101728, 179203,
                                                                       47594, 47720, 106618,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 186628, 0, 3,
                                                                       179203, 101878, 179428,
                                                                       47720, 47846, 106828,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 186943, 0, 3,
                                                                       179428, 102028, 179653,
                                                                       47846, 47972, 107038,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 187258, 0, 3,
                                                                       179653, 102178, 179878,
                                                                       47972, 48098, 107248,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 187573, 0, 3,
                                                                       179878, 102328, 180103,
                                                                       48098, 48224, 107458,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 187888, 0, 3,
                                                                       180328, 102628, 180643,
                                                                       48476, 48644, 107668,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 188308, 0, 3,
                                                                       180643, 102838, 180958,
                                                                       48644, 48812, 107948,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 188728, 0, 3,
                                                                       180958, 103048, 181273,
                                                                       48812, 48980, 108228,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 189148, 0, 3,
                                                                       181273, 103258, 181588,
                                                                       48980, 49148, 108508,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 189568, 0, 3,
                                                                       181588, 103468, 181903,
                                                                       49148, 49316, 108788,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 189988, 0, 3,
                                                                       181903, 103678, 182218,
                                                                       49316, 49484, 109068,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 190408, 0, 3,
                                                                       182218, 103888, 182533,
                                                                       49484, 49652, 109348,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 190828, 0, 3,
                                                                       182533, 104098, 182848,
                                                                       49652, 49820, 109628,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 191248, 0, 3,
                                                                       182848, 104308, 183163,
                                                                       49820, 49988, 109908,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 191668, 0, 3,
                                                                       183163, 104518, 183478,
                                                                       49988, 50156, 110188,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 192088, 0, 3,
                                                                       183478, 104728, 183793,
                                                                       50156, 50324, 110468,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 192508, 0, 3,
                                                                       184108, 105148, 184423,
                                                                       50660, 50828, 110748,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 192928, 0, 3,
                                                                       184423, 105358, 184738,
                                                                       50828, 50996, 111028,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 193348, 0, 3,
                                                                       184738, 105568, 185053,
                                                                       50996, 51164, 111308,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 193768, 0, 3,
                                                                       185053, 105778, 185368,
                                                                       51164, 51332, 111588,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 194188, 0, 3,
                                                                       185368, 105988, 185683,
                                                                       51332, 51500, 111868,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 194608, 0, 3,
                                                                       185683, 106198, 185998,
                                                                       51500, 51668, 112148,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 195028, 0, 3,
                                                                       185998, 106408, 186313,
                                                                       51668, 51836, 112428,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 195448, 0, 3,
                                                                       186313, 106618, 186628,
                                                                       51836, 52004, 112708,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 195868, 0, 3,
                                                                       186628, 106828, 186943,
                                                                       52004, 52172, 112988,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 196288, 0, 3,
                                                                       186943, 107038, 187258,
                                                                       52172, 52340, 113268,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 196708, 0, 3,
                                                                       187258, 107248, 187573,
                                                                       52340, 52508, 113548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 197128, 0, 3,
                                                                       187888, 107668, 188308,
                                                                       52844, 53060, 113828,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 197668, 0, 3,
                                                                       188308, 107948, 188728,
                                                                       53060, 53276, 114188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 198208, 0, 3,
                                                                       188728, 108228, 189148,
                                                                       53276, 53492, 114548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 198748, 0, 3,
                                                                       189148, 108508, 189568,
                                                                       53492, 53708, 114908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 199288, 0, 3,
                                                                       189568, 108788, 189988,
                                                                       53708, 53924, 115268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 199828, 0, 3,
                                                                       189988, 109068, 190408,
                                                                       53924, 54140, 115628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 200368, 0, 3,
                                                                       190408, 109348, 190828,
                                                                       54140, 54356, 115988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 200908, 0, 3,
                                                                       190828, 109628, 191248,
                                                                       54356, 54572, 116348,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 201448, 0, 3,
                                                                       191248, 109908, 191668,
                                                                       54572, 54788, 116708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 201988, 0, 3,
                                                                       191668, 110188, 192088,
                                                                       54788, 55004, 117068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 202528, 0, 3,
                                                                       192508, 110748, 192928,
                                                                       55436, 55652, 117428,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 203068, 0, 3,
                                                                       192928, 111028, 193348,
                                                                       55652, 55868, 117788,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 203608, 0, 3,
                                                                       193348, 111308, 193768,
                                                                       55868, 56084, 118148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 204148, 0, 3,
                                                                       193768, 111588, 194188,
                                                                       56084, 56300, 118508,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 204688, 0, 3,
                                                                       194188, 111868, 194608,
                                                                       56300, 56516, 118868,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 205228, 0, 3,
                                                                       194608, 112148, 195028,
                                                                       56516, 56732, 119228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 205768, 0, 3,
                                                                       195028, 112428, 195448,
                                                                       56732, 56948, 119588,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 206308, 0, 3,
                                                                       195448, 112708, 195868,
                                                                       56948, 57164, 119948,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 206848, 0, 3,
                                                                       195868, 112988, 196288,
                                                                       57164, 57380, 120308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 207388, 0, 3,
                                                                       196288, 113268, 196708,
                                                                       57380, 57596, 120668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 207928, 0, 3,
                                                                       197128, 113828, 197668,
                                                                       58028, 58298, 121028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 208603, 0, 3,
                                                                       197668, 114188, 198208,
                                                                       58298, 58568, 121478,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 209278, 0, 3,
                                                                       198208, 114548, 198748,
                                                                       58568, 58838, 121928,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 209953, 0, 3,
                                                                       198748, 114908, 199288,
                                                                       58838, 59108, 122378,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 210628, 0, 3,
                                                                       199288, 115268, 199828,
                                                                       59108, 59378, 122828,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 211303, 0, 3,
                                                                       199828, 115628, 200368,
                                                                       59378, 59648, 123278,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 211978, 0, 3,
                                                                       200368, 115988, 200908,
                                                                       59648, 59918, 123728,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 212653, 0, 3,
                                                                       200908, 116348, 201448,
                                                                       59918, 60188, 124178,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 213328, 0, 3,
                                                                       201448, 116708, 201988,
                                                                       60188, 60458, 124628,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 214003, 0, 3,
                                                                       202528, 117428, 203068,
                                                                       60998, 61268, 125078,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 214678, 0, 3,
                                                                       203068, 117788, 203608,
                                                                       61268, 61538, 125528,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 215353, 0, 3,
                                                                       203608, 118148, 204148,
                                                                       61538, 61808, 125978,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 216028, 0, 3,
                                                                       204148, 118508, 204688,
                                                                       61808, 62078, 126428,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 216703, 0, 3,
                                                                       204688, 118868, 205228,
                                                                       62078, 62348, 126878,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 217378, 0, 3,
                                                                       205228, 119228, 205768,
                                                                       62348, 62618, 127328,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 218053, 0, 3,
                                                                       205768, 119588, 206308,
                                                                       62618, 62888, 127778,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 218728, 0, 3,
                                                                       206308, 119948, 206848,
                                                                       62888, 63158, 128228,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 219403, 0, 3,
                                                                       206848, 120308, 207388,
                                                                       63158, 63428, 128678,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 220078, 0, 3,
                                                                       207928, 121028, 208603,
                                                                       63968, 64298, 129128,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 220903, 0, 3,
                                                                       208603, 121478, 209278,
                                                                       64298, 64628, 129678,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 221728, 0, 3,
                                                                       209278, 121928, 209953,
                                                                       64628, 64958, 130228,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 222553, 0, 3,
                                                                       209953, 122378, 210628,
                                                                       64958, 65288, 130778,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 223378, 0, 3,
                                                                       210628, 122828, 211303,
                                                                       65288, 65618, 131328,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 224203, 0, 3,
                                                                       211303, 123278, 211978,
                                                                       65618, 65948, 131878,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 225028, 0, 3,
                                                                       211978, 123728, 212653,
                                                                       65948, 66278, 132428,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 225853, 0, 3,
                                                                       212653, 124178, 213328,
                                                                       66278, 66608, 132978,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 226678, 0, 3,
                                                                       214003, 125078, 214678,
                                                                       67268, 67598, 133528,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 227503, 0, 3,
                                                                       214678, 125528, 215353,
                                                                       67598, 67928, 134078,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 228328, 0, 3,
                                                                       215353, 125978, 216028,
                                                                       67928, 68258, 134628,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 229153, 0, 3,
                                                                       216028, 126428, 216703,
                                                                       68258, 68588, 135178,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 229978, 0, 3,
                                                                       216703, 126878, 217378,
                                                                       68588, 68918, 135728,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 230803, 0, 3,
                                                                       217378, 127328, 218053,
                                                                       68918, 69248, 136278,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 231628, 0, 3,
                                                                       218053, 127778, 218728,
                                                                       69248, 69578, 136828,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 232453, 0, 3,
                                                                       218728, 128228, 219403,
                                                                       69578, 69908, 137378,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 233278, 0, 3,
                                                                       220078, 129128, 220903,
                                                                       70568, 70964, 137928,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 234268, 0, 3,
                                                                       220903, 129678, 221728,
                                                                       70964, 71360, 138588,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 235258, 0, 3,
                                                                       221728, 130228, 222553,
                                                                       71360, 71756, 139248,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 236248, 0, 3,
                                                                       222553, 130778, 223378,
                                                                       71756, 72152, 139908,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 237238, 0, 3,
                                                                       223378, 131328, 224203,
                                                                       72152, 72548, 140568,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 238228, 0, 3,
                                                                       224203, 131878, 225028,
                                                                       72548, 72944, 141228,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 239218, 0, 3,
                                                                       225028, 132428, 225853,
                                                                       72944, 73340, 141888,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 240208, 0, 3,
                                                                       226678, 133528, 227503,
                                                                       74132, 74528, 142548,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 241198, 0, 3,
                                                                       227503, 134078, 228328,
                                                                       74528, 74924, 143208,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 242188, 0, 3,
                                                                       228328, 134628, 229153,
                                                                       74924, 75320, 143868,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 243178, 0, 3,
                                                                       229153, 135178, 229978,
                                                                       75320, 75716, 144528,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 244168, 0, 3,
                                                                       229978, 135728, 230803,
                                                                       75716, 76112, 145188,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 245158, 0, 3,
                                                                       230803, 136278, 231628,
                                                                       76112, 76508, 145848,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 246148, 0, 3,
                                                                       231628, 136828, 232453,
                                                                       76508, 76904, 146508,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 247138, 0, 3,
                                                                       233278, 137928, 234268,
                                                                       77696, 78164, 147168,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 248308, 0, 3,
                                                                       234268, 138588, 235258,
                                                                       78164, 78632, 147948,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 249478, 0, 3,
                                                                       235258, 139248, 236248,
                                                                       78632, 79100, 148728,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 250648, 0, 3,
                                                                       236248, 139908, 237238,
                                                                       79100, 79568, 149508,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 251818, 0, 3,
                                                                       237238, 140568, 238228,
                                                                       79568, 80036, 150288,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 252988, 0, 3,
                                                                       238228, 141228, 239218,
                                                                       80036, 80504, 151068,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 254158, 0, 3,
                                                                       240208, 142548, 241198,
                                                                       81440, 81908, 151848,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 255328, 0, 3,
                                                                       241198, 143208, 242188,
                                                                       81908, 82376, 152628,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 256498, 0, 3,
                                                                       242188, 143868, 243178,
                                                                       82376, 82844, 153408,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 257668, 0, 3,
                                                                       243178, 144528, 244168,
                                                                       82844, 83312, 154188,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 258838, 0, 3,
                                                                       244168, 145188, 245158,
                                                                       83312, 83780, 154968,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 260008, 0, 3,
                                                                       245158, 145848, 246148,
                                                                       83780, 84248, 155748,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 261178, 0, 3,
                                                                       247138, 147168, 248308,
                                                                       85184, 85730, 156528,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 262543, 0, 3,
                                                                       248308, 147948, 249478,
                                                                       85730, 86276, 157438,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 263908, 0, 3,
                                                                       249478, 148728, 250648,
                                                                       86276, 86822, 158348,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 265273, 0, 3,
                                                                       250648, 149508, 251818,
                                                                       86822, 87368, 159258,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 266638, 0, 3,
                                                                       251818, 150288, 252988,
                                                                       87368, 87914, 160168,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 268003, 0, 3,
                                                                       254158, 151848, 255328,
                                                                       89006, 89552, 161078,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 269368, 0, 3,
                                                                       255328, 152628, 256498,
                                                                       89552, 90098, 161988,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 270733, 0, 3,
                                                                       256498, 153408, 257668,
                                                                       90098, 90644, 162898,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 272098, 0, 3,
                                                                       257668, 154188, 258838,
                                                                       90644, 91190, 163808,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 273463, 0, 3,
                                                                       258838, 154968, 260008,
                                                                       91190, 91736, 164718,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 274828, 3, 92828,
                                                                       92838, 165658, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 274849, 3, 92838,
                                                                       92848, 165673, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 274870, 3, 92848,
                                                                       92858, 165688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 274891, 3, 92858,
                                                                       92868, 165703, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 274912, 3, 92868,
                                                                       92878, 165718, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 274933, 3, 92878,
                                                                       92888, 165733, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 274954, 3, 92888,
                                                                       92898, 165748, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 274975, 3, 92898,
                                                                       92908, 165763, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 274996, 3, 92908,
                                                                       92918, 165778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275017, 3, 92918,
                                                                       92928, 165793, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275038, 3, 92928,
                                                                       92938, 165808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275059, 3, 92938,
                                                                       92948, 165823, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275080, 3, 92948,
                                                                       92958, 165838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275101, 3, 92958,
                                                                       92968, 165853, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275122, 3, 92968,
                                                                       92978, 165868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275143, 3, 92998,
                                                                       93008, 165913, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275164, 3, 93008,
                                                                       93018, 165928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275185, 3, 93018,
                                                                       93028, 165943, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275206, 3, 93028,
                                                                       93038, 165958, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275227, 3, 93038,
                                                                       93048, 165973, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275248, 3, 93048,
                                                                       93058, 165988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275269, 3, 93058,
                                                                       93068, 166003, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275290, 3, 93068,
                                                                       93078, 166018, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275311, 3, 93078,
                                                                       93088, 166033, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275332, 3, 93088,
                                                                       93098, 166048, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275353, 3, 93098,
                                                                       93108, 166063, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275374, 3, 93108,
                                                                       93118, 166078, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275395, 3, 93118,
                                                                       93128, 166093, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275416, 3, 93128,
                                                                       93138, 166108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 275437, 3, 93138,
                                                                       93148, 166123, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 275458, 0, 3,
                                                                       274828, 165658, 274849,
                                                                       93168, 93198, 166228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 275521, 0, 3,
                                                                       274849, 165673, 274870,
                                                                       93198, 93228, 166273,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 275584, 0, 3,
                                                                       274870, 165688, 274891,
                                                                       93228, 93258, 166318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 275647, 0, 3,
                                                                       274891, 165703, 274912,
                                                                       93258, 93288, 166363,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 275710, 0, 3,
                                                                       274912, 165718, 274933,
                                                                       93288, 93318, 166408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 275773, 0, 3,
                                                                       274933, 165733, 274954,
                                                                       93318, 93348, 166453,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 275836, 0, 3,
                                                                       274954, 165748, 274975,
                                                                       93348, 93378, 166498,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 275899, 0, 3,
                                                                       274975, 165763, 274996,
                                                                       93378, 93408, 166543,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 275962, 0, 3,
                                                                       274996, 165778, 275017,
                                                                       93408, 93438, 166588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276025, 0, 3,
                                                                       275017, 165793, 275038,
                                                                       93438, 93468, 166633,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276088, 0, 3,
                                                                       275038, 165808, 275059,
                                                                       93468, 93498, 166678,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276151, 0, 3,
                                                                       275059, 165823, 275080,
                                                                       93498, 93528, 166723,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276214, 0, 3,
                                                                       275080, 165838, 275101,
                                                                       93528, 93558, 166768,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276277, 0, 3,
                                                                       275101, 165853, 275122,
                                                                       93558, 93588, 166813,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276340, 0, 3,
                                                                       275143, 165913, 275164,
                                                                       93648, 93678, 166948,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276403, 0, 3,
                                                                       275164, 165928, 275185,
                                                                       93678, 93708, 166993,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276466, 0, 3,
                                                                       275185, 165943, 275206,
                                                                       93708, 93738, 167038,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276529, 0, 3,
                                                                       275206, 165958, 275227,
                                                                       93738, 93768, 167083,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276592, 0, 3,
                                                                       275227, 165973, 275248,
                                                                       93768, 93798, 167128,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276655, 0, 3,
                                                                       275248, 165988, 275269,
                                                                       93798, 93828, 167173,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276718, 0, 3,
                                                                       275269, 166003, 275290,
                                                                       93828, 93858, 167218,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276781, 0, 3,
                                                                       275290, 166018, 275311,
                                                                       93858, 93888, 167263,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276844, 0, 3,
                                                                       275311, 166033, 275332,
                                                                       93888, 93918, 167308,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276907, 0, 3,
                                                                       275332, 166048, 275353,
                                                                       93918, 93948, 167353,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 276970, 0, 3,
                                                                       275353, 166063, 275374,
                                                                       93948, 93978, 167398,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 277033, 0, 3,
                                                                       275374, 166078, 275395,
                                                                       93978, 94008, 167443,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 277096, 0, 3,
                                                                       275395, 166093, 275416,
                                                                       94008, 94038, 167488,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 277159, 0, 3,
                                                                       275416, 166108, 275437,
                                                                       94038, 94068, 167533,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 277222, 0, 3,
                                                                       275458, 166228, 275521,
                                                                       94128, 94188, 167758,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 277348, 0, 3,
                                                                       275521, 166273, 275584,
                                                                       94188, 94248, 167848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 277474, 0, 3,
                                                                       275584, 166318, 275647,
                                                                       94248, 94308, 167938,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 277600, 0, 3,
                                                                       275647, 166363, 275710,
                                                                       94308, 94368, 168028,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 277726, 0, 3,
                                                                       275710, 166408, 275773,
                                                                       94368, 94428, 168118,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 277852, 0, 3,
                                                                       275773, 166453, 275836,
                                                                       94428, 94488, 168208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 277978, 0, 3,
                                                                       275836, 166498, 275899,
                                                                       94488, 94548, 168298,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 278104, 0, 3,
                                                                       275899, 166543, 275962,
                                                                       94548, 94608, 168388,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 278230, 0, 3,
                                                                       275962, 166588, 276025,
                                                                       94608, 94668, 168478,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 278356, 0, 3,
                                                                       276025, 166633, 276088,
                                                                       94668, 94728, 168568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 278482, 0, 3,
                                                                       276088, 166678, 276151,
                                                                       94728, 94788, 168658,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 278608, 0, 3,
                                                                       276151, 166723, 276214,
                                                                       94788, 94848, 168748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 278734, 0, 3,
                                                                       276214, 166768, 276277,
                                                                       94848, 94908, 168838,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 278860, 0, 3,
                                                                       276340, 166948, 276403,
                                                                       95028, 95088, 169108,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 278986, 0, 3,
                                                                       276403, 166993, 276466,
                                                                       95088, 95148, 169198,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 279112, 0, 3,
                                                                       276466, 167038, 276529,
                                                                       95148, 95208, 169288,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 279238, 0, 3,
                                                                       276529, 167083, 276592,
                                                                       95208, 95268, 169378,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 279364, 0, 3,
                                                                       276592, 167128, 276655,
                                                                       95268, 95328, 169468,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 279490, 0, 3,
                                                                       276655, 167173, 276718,
                                                                       95328, 95388, 169558,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 279616, 0, 3,
                                                                       276718, 167218, 276781,
                                                                       95388, 95448, 169648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 279742, 0, 3,
                                                                       276781, 167263, 276844,
                                                                       95448, 95508, 169738,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 279868, 0, 3,
                                                                       276844, 167308, 276907,
                                                                       95508, 95568, 169828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 279994, 0, 3,
                                                                       276907, 167353, 276970,
                                                                       95568, 95628, 169918,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 280120, 0, 3,
                                                                       276970, 167398, 277033,
                                                                       95628, 95688, 170008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 280246, 0, 3,
                                                                       277033, 167443, 277096,
                                                                       95688, 95748, 170098,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 280372, 0, 3,
                                                                       277096, 167488, 277159,
                                                                       95748, 95808, 170188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 280498, 0, 3,
                                                                       277222, 167758, 277348,
                                                                       95928, 96028, 170578,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 280708, 0, 3,
                                                                       277348, 167848, 277474,
                                                                       96028, 96128, 170728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 280918, 0, 3,
                                                                       277474, 167938, 277600,
                                                                       96128, 96228, 170878,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 281128, 0, 3,
                                                                       277600, 168028, 277726,
                                                                       96228, 96328, 171028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 281338, 0, 3,
                                                                       277726, 168118, 277852,
                                                                       96328, 96428, 171178,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 281548, 0, 3,
                                                                       277852, 168208, 277978,
                                                                       96428, 96528, 171328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 281758, 0, 3,
                                                                       277978, 168298, 278104,
                                                                       96528, 96628, 171478,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 281968, 0, 3,
                                                                       278104, 168388, 278230,
                                                                       96628, 96728, 171628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 282178, 0, 3,
                                                                       278230, 168478, 278356,
                                                                       96728, 96828, 171778,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 282388, 0, 3,
                                                                       278356, 168568, 278482,
                                                                       96828, 96928, 171928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 282598, 0, 3,
                                                                       278482, 168658, 278608,
                                                                       96928, 97028, 172078,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 282808, 0, 3,
                                                                       278608, 168748, 278734,
                                                                       97028, 97128, 172228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 283018, 0, 3,
                                                                       278860, 169108, 278986,
                                                                       97328, 97428, 172678,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 283228, 0, 3,
                                                                       278986, 169198, 279112,
                                                                       97428, 97528, 172828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 283438, 0, 3,
                                                                       279112, 169288, 279238,
                                                                       97528, 97628, 172978,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 283648, 0, 3,
                                                                       279238, 169378, 279364,
                                                                       97628, 97728, 173128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 283858, 0, 3,
                                                                       279364, 169468, 279490,
                                                                       97728, 97828, 173278,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 284068, 0, 3,
                                                                       279490, 169558, 279616,
                                                                       97828, 97928, 173428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 284278, 0, 3,
                                                                       279616, 169648, 279742,
                                                                       97928, 98028, 173578,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 284488, 0, 3,
                                                                       279742, 169738, 279868,
                                                                       98028, 98128, 173728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 284698, 0, 3,
                                                                       279868, 169828, 279994,
                                                                       98128, 98228, 173878,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 284908, 0, 3,
                                                                       279994, 169918, 280120,
                                                                       98228, 98328, 174028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 285118, 0, 3,
                                                                       280120, 170008, 280246,
                                                                       98328, 98428, 174178,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 285328, 0, 3,
                                                                       280246, 170098, 280372,
                                                                       98428, 98528, 174328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 285538, 0, 3,
                                                                       280498, 170578, 280708,
                                                                       98728, 98878, 174928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 285853, 0, 3,
                                                                       280708, 170728, 280918,
                                                                       98878, 99028, 175153,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 286168, 0, 3,
                                                                       280918, 170878, 281128,
                                                                       99028, 99178, 175378,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 286483, 0, 3,
                                                                       281128, 171028, 281338,
                                                                       99178, 99328, 175603,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 286798, 0, 3,
                                                                       281338, 171178, 281548,
                                                                       99328, 99478, 175828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 287113, 0, 3,
                                                                       281548, 171328, 281758,
                                                                       99478, 99628, 176053,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 287428, 0, 3,
                                                                       281758, 171478, 281968,
                                                                       99628, 99778, 176278,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 287743, 0, 3,
                                                                       281968, 171628, 282178,
                                                                       99778, 99928, 176503,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 288058, 0, 3,
                                                                       282178, 171778, 282388,
                                                                       99928, 100078, 176728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 288373, 0, 3,
                                                                       282388, 171928, 282598,
                                                                       100078, 100228, 176953,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 288688, 0, 3,
                                                                       282598, 172078, 282808,
                                                                       100228, 100378, 177178,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 289003, 0, 3,
                                                                       283018, 172678, 283228,
                                                                       100678, 100828, 177853,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 289318, 0, 3,
                                                                       283228, 172828, 283438,
                                                                       100828, 100978, 178078,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 289633, 0, 3,
                                                                       283438, 172978, 283648,
                                                                       100978, 101128, 178303,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 289948, 0, 3,
                                                                       283648, 173128, 283858,
                                                                       101128, 101278, 178528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 290263, 0, 3,
                                                                       283858, 173278, 284068,
                                                                       101278, 101428, 178753,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 290578, 0, 3,
                                                                       284068, 173428, 284278,
                                                                       101428, 101578, 178978,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 290893, 0, 3,
                                                                       284278, 173578, 284488,
                                                                       101578, 101728, 179203,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 291208, 0, 3,
                                                                       284488, 173728, 284698,
                                                                       101728, 101878, 179428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 291523, 0, 3,
                                                                       284698, 173878, 284908,
                                                                       101878, 102028, 179653,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 291838, 0, 3,
                                                                       284908, 174028, 285118,
                                                                       102028, 102178, 179878,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 292153, 0, 3,
                                                                       285118, 174178, 285328,
                                                                       102178, 102328, 180103,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 292468, 0, 3,
                                                                       285538, 174928, 285853,
                                                                       102628, 102838, 180958,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 292909, 0, 3,
                                                                       285853, 175153, 286168,
                                                                       102838, 103048, 181273,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 293350, 0, 3,
                                                                       286168, 175378, 286483,
                                                                       103048, 103258, 181588,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 293791, 0, 3,
                                                                       286483, 175603, 286798,
                                                                       103258, 103468, 181903,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 294232, 0, 3,
                                                                       286798, 175828, 287113,
                                                                       103468, 103678, 182218,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 294673, 0, 3,
                                                                       287113, 176053, 287428,
                                                                       103678, 103888, 182533,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 295114, 0, 3,
                                                                       287428, 176278, 287743,
                                                                       103888, 104098, 182848,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 295555, 0, 3,
                                                                       287743, 176503, 288058,
                                                                       104098, 104308, 183163,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 295996, 0, 3,
                                                                       288058, 176728, 288373,
                                                                       104308, 104518, 183478,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 296437, 0, 3,
                                                                       288373, 176953, 288688,
                                                                       104518, 104728, 183793,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 296878, 0, 3,
                                                                       289003, 177853, 289318,
                                                                       105148, 105358, 184738,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 297319, 0, 3,
                                                                       289318, 178078, 289633,
                                                                       105358, 105568, 185053,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 297760, 0, 3,
                                                                       289633, 178303, 289948,
                                                                       105568, 105778, 185368,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 298201, 0, 3,
                                                                       289948, 178528, 290263,
                                                                       105778, 105988, 185683,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 298642, 0, 3,
                                                                       290263, 178753, 290578,
                                                                       105988, 106198, 185998,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 299083, 0, 3,
                                                                       290578, 178978, 290893,
                                                                       106198, 106408, 186313,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 299524, 0, 3,
                                                                       290893, 179203, 291208,
                                                                       106408, 106618, 186628,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 299965, 0, 3,
                                                                       291208, 179428, 291523,
                                                                       106618, 106828, 186943,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 300406, 0, 3,
                                                                       291523, 179653, 291838,
                                                                       106828, 107038, 187258,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 300847, 0, 3,
                                                                       291838, 179878, 292153,
                                                                       107038, 107248, 187573,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 301288, 0, 3,
                                                                       292468, 180958, 292909,
                                                                       107668, 107948, 188728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 301876, 0, 3,
                                                                       292909, 181273, 293350,
                                                                       107948, 108228, 189148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 302464, 0, 3,
                                                                       293350, 181588, 293791,
                                                                       108228, 108508, 189568,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 303052, 0, 3,
                                                                       293791, 181903, 294232,
                                                                       108508, 108788, 189988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 303640, 0, 3,
                                                                       294232, 182218, 294673,
                                                                       108788, 109068, 190408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 304228, 0, 3,
                                                                       294673, 182533, 295114,
                                                                       109068, 109348, 190828,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 304816, 0, 3,
                                                                       295114, 182848, 295555,
                                                                       109348, 109628, 191248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 305404, 0, 3,
                                                                       295555, 183163, 295996,
                                                                       109628, 109908, 191668,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 305992, 0, 3,
                                                                       295996, 183478, 296437,
                                                                       109908, 110188, 192088,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 306580, 0, 3,
                                                                       296878, 184738, 297319,
                                                                       110748, 111028, 193348,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 307168, 0, 3,
                                                                       297319, 185053, 297760,
                                                                       111028, 111308, 193768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 307756, 0, 3,
                                                                       297760, 185368, 298201,
                                                                       111308, 111588, 194188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 308344, 0, 3,
                                                                       298201, 185683, 298642,
                                                                       111588, 111868, 194608,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 308932, 0, 3,
                                                                       298642, 185998, 299083,
                                                                       111868, 112148, 195028,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 309520, 0, 3,
                                                                       299083, 186313, 299524,
                                                                       112148, 112428, 195448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 310108, 0, 3,
                                                                       299524, 186628, 299965,
                                                                       112428, 112708, 195868,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 310696, 0, 3,
                                                                       299965, 186943, 300406,
                                                                       112708, 112988, 196288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 311284, 0, 3,
                                                                       300406, 187258, 300847,
                                                                       112988, 113268, 196708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 311872, 0, 3,
                                                                       301288, 188728, 301876,
                                                                       113828, 114188, 198208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 312628, 0, 3,
                                                                       301876, 189148, 302464,
                                                                       114188, 114548, 198748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 313384, 0, 3,
                                                                       302464, 189568, 303052,
                                                                       114548, 114908, 199288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 314140, 0, 3,
                                                                       303052, 189988, 303640,
                                                                       114908, 115268, 199828,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 314896, 0, 3,
                                                                       303640, 190408, 304228,
                                                                       115268, 115628, 200368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 315652, 0, 3,
                                                                       304228, 190828, 304816,
                                                                       115628, 115988, 200908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 316408, 0, 3,
                                                                       304816, 191248, 305404,
                                                                       115988, 116348, 201448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 317164, 0, 3,
                                                                       305404, 191668, 305992,
                                                                       116348, 116708, 201988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 317920, 0, 3,
                                                                       306580, 193348, 307168,
                                                                       117428, 117788, 203608,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 318676, 0, 3,
                                                                       307168, 193768, 307756,
                                                                       117788, 118148, 204148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 319432, 0, 3,
                                                                       307756, 194188, 308344,
                                                                       118148, 118508, 204688,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 320188, 0, 3,
                                                                       308344, 194608, 308932,
                                                                       118508, 118868, 205228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 320944, 0, 3,
                                                                       308932, 195028, 309520,
                                                                       118868, 119228, 205768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 321700, 0, 3,
                                                                       309520, 195448, 310108,
                                                                       119228, 119588, 206308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 322456, 0, 3,
                                                                       310108, 195868, 310696,
                                                                       119588, 119948, 206848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 323212, 0, 3,
                                                                       310696, 196288, 311284,
                                                                       119948, 120308, 207388,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 323968, 0, 3,
                                                                       311872, 198208, 312628,
                                                                       121028, 121478, 209278,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 324913, 0, 3,
                                                                       312628, 198748, 313384,
                                                                       121478, 121928, 209953,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 325858, 0, 3,
                                                                       313384, 199288, 314140,
                                                                       121928, 122378, 210628,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 326803, 0, 3,
                                                                       314140, 199828, 314896,
                                                                       122378, 122828, 211303,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 327748, 0, 3,
                                                                       314896, 200368, 315652,
                                                                       122828, 123278, 211978,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 328693, 0, 3,
                                                                       315652, 200908, 316408,
                                                                       123278, 123728, 212653,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 329638, 0, 3,
                                                                       316408, 201448, 317164,
                                                                       123728, 124178, 213328,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 330583, 0, 3,
                                                                       317920, 203608, 318676,
                                                                       125078, 125528, 215353,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 331528, 0, 3,
                                                                       318676, 204148, 319432,
                                                                       125528, 125978, 216028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 332473, 0, 3,
                                                                       319432, 204688, 320188,
                                                                       125978, 126428, 216703,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 333418, 0, 3,
                                                                       320188, 205228, 320944,
                                                                       126428, 126878, 217378,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 334363, 0, 3,
                                                                       320944, 205768, 321700,
                                                                       126878, 127328, 218053,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 335308, 0, 3,
                                                                       321700, 206308, 322456,
                                                                       127328, 127778, 218728,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 336253, 0, 3,
                                                                       322456, 206848, 323212,
                                                                       127778, 128228, 219403,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 337198, 0, 3,
                                                                       323968, 209278, 324913,
                                                                       129128, 129678, 221728,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 338353, 0, 3,
                                                                       324913, 209953, 325858,
                                                                       129678, 130228, 222553,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 339508, 0, 3,
                                                                       325858, 210628, 326803,
                                                                       130228, 130778, 223378,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 340663, 0, 3,
                                                                       326803, 211303, 327748,
                                                                       130778, 131328, 224203,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 341818, 0, 3,
                                                                       327748, 211978, 328693,
                                                                       131328, 131878, 225028,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 342973, 0, 3,
                                                                       328693, 212653, 329638,
                                                                       131878, 132428, 225853,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 344128, 0, 3,
                                                                       330583, 215353, 331528,
                                                                       133528, 134078, 228328,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 345283, 0, 3,
                                                                       331528, 216028, 332473,
                                                                       134078, 134628, 229153,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 346438, 0, 3,
                                                                       332473, 216703, 333418,
                                                                       134628, 135178, 229978,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 347593, 0, 3,
                                                                       333418, 217378, 334363,
                                                                       135178, 135728, 230803,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 348748, 0, 3,
                                                                       334363, 218053, 335308,
                                                                       135728, 136278, 231628,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 349903, 0, 3,
                                                                       335308, 218728, 336253,
                                                                       136278, 136828, 232453,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 351058, 0, 3,
                                                                       337198, 221728, 338353,
                                                                       137928, 138588, 235258,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 352444, 0, 3,
                                                                       338353, 222553, 339508,
                                                                       138588, 139248, 236248,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 353830, 0, 3,
                                                                       339508, 223378, 340663,
                                                                       139248, 139908, 237238,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 355216, 0, 3,
                                                                       340663, 224203, 341818,
                                                                       139908, 140568, 238228,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 356602, 0, 3,
                                                                       341818, 225028, 342973,
                                                                       140568, 141228, 239218,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 357988, 0, 3,
                                                                       344128, 228328, 345283,
                                                                       142548, 143208, 242188,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 359374, 0, 3,
                                                                       345283, 229153, 346438,
                                                                       143208, 143868, 243178,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 360760, 0, 3,
                                                                       346438, 229978, 347593,
                                                                       143868, 144528, 244168,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 362146, 0, 3,
                                                                       347593, 230803, 348748,
                                                                       144528, 145188, 245158,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 363532, 0, 3,
                                                                       348748, 231628, 349903,
                                                                       145188, 145848, 246148,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 364918, 0, 3,
                                                                       351058, 235258, 352444,
                                                                       147168, 147948, 249478,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 366556, 0, 3,
                                                                       352444, 236248, 353830,
                                                                       147948, 148728, 250648,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 368194, 0, 3,
                                                                       353830, 237238, 355216,
                                                                       148728, 149508, 251818,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 369832, 0, 3,
                                                                       355216, 238228, 356602,
                                                                       149508, 150288, 252988,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 371470, 0, 3,
                                                                       357988, 242188, 359374,
                                                                       151848, 152628, 256498,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 373108, 0, 3,
                                                                       359374, 243178, 360760,
                                                                       152628, 153408, 257668,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 374746, 0, 3,
                                                                       360760, 244168, 362146,
                                                                       153408, 154188, 258838,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 376384, 0, 3,
                                                                       362146, 245158, 363532,
                                                                       154188, 154968, 260008,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 378022, 0, 3,
                                                                       364918, 249478, 366556,
                                                                       156528, 157438, 263908,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 379933, 0, 3,
                                                                       366556, 250648, 368194,
                                                                       157438, 158348, 265273,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 381844, 0, 3,
                                                                       368194, 251818, 369832,
                                                                       158348, 159258, 266638,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 383755, 0, 3,
                                                                       371470, 256498, 373108,
                                                                       161078, 161988, 270733,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 385666, 0, 3,
                                                                       373108, 257668, 374746,
                                                                       161988, 162898, 272098,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 387577, 0, 3,
                                                                       374746, 258838, 376384,
                                                                       162898, 163808, 273463,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389488, 3, 165628,
                                                                       165643, 274828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389516, 3, 165643,
                                                                       165658, 274849, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389544, 3, 165658,
                                                                       165673, 274870, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389572, 3, 165673,
                                                                       165688, 274891, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389600, 3, 165688,
                                                                       165703, 274912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389628, 3, 165703,
                                                                       165718, 274933, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389656, 3, 165718,
                                                                       165733, 274954, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389684, 3, 165733,
                                                                       165748, 274975, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389712, 3, 165748,
                                                                       165763, 274996, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389740, 3, 165763,
                                                                       165778, 275017, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389768, 3, 165778,
                                                                       165793, 275038, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389796, 3, 165793,
                                                                       165808, 275059, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389824, 3, 165808,
                                                                       165823, 275080, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389852, 3, 165823,
                                                                       165838, 275101, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389880, 3, 165838,
                                                                       165853, 275122, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389908, 3, 165883,
                                                                       165898, 275143, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389936, 3, 165898,
                                                                       165913, 275164, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389964, 3, 165913,
                                                                       165928, 275185, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 389992, 3, 165928,
                                                                       165943, 275206, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 390020, 3, 165943,
                                                                       165958, 275227, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 390048, 3, 165958,
                                                                       165973, 275248, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 390076, 3, 165973,
                                                                       165988, 275269, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 390104, 3, 165988,
                                                                       166003, 275290, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 390132, 3, 166003,
                                                                       166018, 275311, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 390160, 3, 166018,
                                                                       166033, 275332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 390188, 3, 166033,
                                                                       166048, 275353, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 390216, 3, 166048,
                                                                       166063, 275374, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 390244, 3, 166063,
                                                                       166078, 275395, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 390272, 3, 166078,
                                                                       166093, 275416, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 390300, 3, 166093,
                                                                       166108, 275437, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 390328, 0, 3,
                                                                       389488, 274828, 389516,
                                                                       166138, 166183, 275458,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 390412, 0, 3,
                                                                       389516, 274849, 389544,
                                                                       166183, 166228, 275521,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 390496, 0, 3,
                                                                       389544, 274870, 389572,
                                                                       166228, 166273, 275584,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 390580, 0, 3,
                                                                       389572, 274891, 389600,
                                                                       166273, 166318, 275647,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 390664, 0, 3,
                                                                       389600, 274912, 389628,
                                                                       166318, 166363, 275710,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 390748, 0, 3,
                                                                       389628, 274933, 389656,
                                                                       166363, 166408, 275773,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 390832, 0, 3,
                                                                       389656, 274954, 389684,
                                                                       166408, 166453, 275836,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 390916, 0, 3,
                                                                       389684, 274975, 389712,
                                                                       166453, 166498, 275899,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391000, 0, 3,
                                                                       389712, 274996, 389740,
                                                                       166498, 166543, 275962,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391084, 0, 3,
                                                                       389740, 275017, 389768,
                                                                       166543, 166588, 276025,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391168, 0, 3,
                                                                       389768, 275038, 389796,
                                                                       166588, 166633, 276088,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391252, 0, 3,
                                                                       389796, 275059, 389824,
                                                                       166633, 166678, 276151,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391336, 0, 3,
                                                                       389824, 275080, 389852,
                                                                       166678, 166723, 276214,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391420, 0, 3,
                                                                       389852, 275101, 389880,
                                                                       166723, 166768, 276277,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391504, 0, 3,
                                                                       389908, 275143, 389936,
                                                                       166858, 166903, 276340,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391588, 0, 3,
                                                                       389936, 275164, 389964,
                                                                       166903, 166948, 276403,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391672, 0, 3,
                                                                       389964, 275185, 389992,
                                                                       166948, 166993, 276466,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391756, 0, 3,
                                                                       389992, 275206, 390020,
                                                                       166993, 167038, 276529,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391840, 0, 3,
                                                                       390020, 275227, 390048,
                                                                       167038, 167083, 276592,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 391924, 0, 3,
                                                                       390048, 275248, 390076,
                                                                       167083, 167128, 276655,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 392008, 0, 3,
                                                                       390076, 275269, 390104,
                                                                       167128, 167173, 276718,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 392092, 0, 3,
                                                                       390104, 275290, 390132,
                                                                       167173, 167218, 276781,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 392176, 0, 3,
                                                                       390132, 275311, 390160,
                                                                       167218, 167263, 276844,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 392260, 0, 3,
                                                                       390160, 275332, 390188,
                                                                       167263, 167308, 276907,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 392344, 0, 3,
                                                                       390188, 275353, 390216,
                                                                       167308, 167353, 276970,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 392428, 0, 3,
                                                                       390216, 275374, 390244,
                                                                       167353, 167398, 277033,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 392512, 0, 3,
                                                                       390244, 275395, 390272,
                                                                       167398, 167443, 277096,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 392596, 0, 3,
                                                                       390272, 275416, 390300,
                                                                       167443, 167488, 277159,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 392680, 0, 3,
                                                                       390328, 275458, 390412,
                                                                       167578, 167668, 277222,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 392848, 0, 3,
                                                                       390412, 275521, 390496,
                                                                       167668, 167758, 277348,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 393016, 0, 3,
                                                                       390496, 275584, 390580,
                                                                       167758, 167848, 277474,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 393184, 0, 3,
                                                                       390580, 275647, 390664,
                                                                       167848, 167938, 277600,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 393352, 0, 3,
                                                                       390664, 275710, 390748,
                                                                       167938, 168028, 277726,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 393520, 0, 3,
                                                                       390748, 275773, 390832,
                                                                       168028, 168118, 277852,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 393688, 0, 3,
                                                                       390832, 275836, 390916,
                                                                       168118, 168208, 277978,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 393856, 0, 3,
                                                                       390916, 275899, 391000,
                                                                       168208, 168298, 278104,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 394024, 0, 3,
                                                                       391000, 275962, 391084,
                                                                       168298, 168388, 278230,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 394192, 0, 3,
                                                                       391084, 276025, 391168,
                                                                       168388, 168478, 278356,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 394360, 0, 3,
                                                                       391168, 276088, 391252,
                                                                       168478, 168568, 278482,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 394528, 0, 3,
                                                                       391252, 276151, 391336,
                                                                       168568, 168658, 278608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 394696, 0, 3,
                                                                       391336, 276214, 391420,
                                                                       168658, 168748, 278734,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 394864, 0, 3,
                                                                       391504, 276340, 391588,
                                                                       168928, 169018, 278860,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 395032, 0, 3,
                                                                       391588, 276403, 391672,
                                                                       169018, 169108, 278986,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 395200, 0, 3,
                                                                       391672, 276466, 391756,
                                                                       169108, 169198, 279112,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 395368, 0, 3,
                                                                       391756, 276529, 391840,
                                                                       169198, 169288, 279238,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 395536, 0, 3,
                                                                       391840, 276592, 391924,
                                                                       169288, 169378, 279364,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 395704, 0, 3,
                                                                       391924, 276655, 392008,
                                                                       169378, 169468, 279490,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 395872, 0, 3,
                                                                       392008, 276718, 392092,
                                                                       169468, 169558, 279616,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 396040, 0, 3,
                                                                       392092, 276781, 392176,
                                                                       169558, 169648, 279742,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 396208, 0, 3,
                                                                       392176, 276844, 392260,
                                                                       169648, 169738, 279868,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 396376, 0, 3,
                                                                       392260, 276907, 392344,
                                                                       169738, 169828, 279994,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 396544, 0, 3,
                                                                       392344, 276970, 392428,
                                                                       169828, 169918, 280120,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 396712, 0, 3,
                                                                       392428, 277033, 392512,
                                                                       169918, 170008, 280246,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 396880, 0, 3,
                                                                       392512, 277096, 392596,
                                                                       170008, 170098, 280372,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 397048, 0, 3,
                                                                       392680, 277222, 392848,
                                                                       170278, 170428, 280498,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 397328, 0, 3,
                                                                       392848, 277348, 393016,
                                                                       170428, 170578, 280708,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 397608, 0, 3,
                                                                       393016, 277474, 393184,
                                                                       170578, 170728, 280918,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 397888, 0, 3,
                                                                       393184, 277600, 393352,
                                                                       170728, 170878, 281128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 398168, 0, 3,
                                                                       393352, 277726, 393520,
                                                                       170878, 171028, 281338,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 398448, 0, 3,
                                                                       393520, 277852, 393688,
                                                                       171028, 171178, 281548,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 398728, 0, 3,
                                                                       393688, 277978, 393856,
                                                                       171178, 171328, 281758,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 399008, 0, 3,
                                                                       393856, 278104, 394024,
                                                                       171328, 171478, 281968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 399288, 0, 3,
                                                                       394024, 278230, 394192,
                                                                       171478, 171628, 282178,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 399568, 0, 3,
                                                                       394192, 278356, 394360,
                                                                       171628, 171778, 282388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 399848, 0, 3,
                                                                       394360, 278482, 394528,
                                                                       171778, 171928, 282598,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 400128, 0, 3,
                                                                       394528, 278608, 394696,
                                                                       171928, 172078, 282808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 400408, 0, 3,
                                                                       394864, 278860, 395032,
                                                                       172378, 172528, 283018,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 400688, 0, 3,
                                                                       395032, 278986, 395200,
                                                                       172528, 172678, 283228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 400968, 0, 3,
                                                                       395200, 279112, 395368,
                                                                       172678, 172828, 283438,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 401248, 0, 3,
                                                                       395368, 279238, 395536,
                                                                       172828, 172978, 283648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 401528, 0, 3,
                                                                       395536, 279364, 395704,
                                                                       172978, 173128, 283858,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 401808, 0, 3,
                                                                       395704, 279490, 395872,
                                                                       173128, 173278, 284068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 402088, 0, 3,
                                                                       395872, 279616, 396040,
                                                                       173278, 173428, 284278,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 402368, 0, 3,
                                                                       396040, 279742, 396208,
                                                                       173428, 173578, 284488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 402648, 0, 3,
                                                                       396208, 279868, 396376,
                                                                       173578, 173728, 284698,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 402928, 0, 3,
                                                                       396376, 279994, 396544,
                                                                       173728, 173878, 284908,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 403208, 0, 3,
                                                                       396544, 280120, 396712,
                                                                       173878, 174028, 285118,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 403488, 0, 3,
                                                                       396712, 280246, 396880,
                                                                       174028, 174178, 285328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 403768, 0, 3,
                                                                       397048, 280498, 397328,
                                                                       174478, 174703, 285538,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 404188, 0, 3,
                                                                       397328, 280708, 397608,
                                                                       174703, 174928, 285853,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 404608, 0, 3,
                                                                       397608, 280918, 397888,
                                                                       174928, 175153, 286168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 405028, 0, 3,
                                                                       397888, 281128, 398168,
                                                                       175153, 175378, 286483,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 405448, 0, 3,
                                                                       398168, 281338, 398448,
                                                                       175378, 175603, 286798,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 405868, 0, 3,
                                                                       398448, 281548, 398728,
                                                                       175603, 175828, 287113,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 406288, 0, 3,
                                                                       398728, 281758, 399008,
                                                                       175828, 176053, 287428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 406708, 0, 3,
                                                                       399008, 281968, 399288,
                                                                       176053, 176278, 287743,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 407128, 0, 3,
                                                                       399288, 282178, 399568,
                                                                       176278, 176503, 288058,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 407548, 0, 3,
                                                                       399568, 282388, 399848,
                                                                       176503, 176728, 288373,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 407968, 0, 3,
                                                                       399848, 282598, 400128,
                                                                       176728, 176953, 288688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 408388, 0, 3,
                                                                       400408, 283018, 400688,
                                                                       177403, 177628, 289003,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 408808, 0, 3,
                                                                       400688, 283228, 400968,
                                                                       177628, 177853, 289318,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 409228, 0, 3,
                                                                       400968, 283438, 401248,
                                                                       177853, 178078, 289633,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 409648, 0, 3,
                                                                       401248, 283648, 401528,
                                                                       178078, 178303, 289948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 410068, 0, 3,
                                                                       401528, 283858, 401808,
                                                                       178303, 178528, 290263,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 410488, 0, 3,
                                                                       401808, 284068, 402088,
                                                                       178528, 178753, 290578,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 410908, 0, 3,
                                                                       402088, 284278, 402368,
                                                                       178753, 178978, 290893,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 411328, 0, 3,
                                                                       402368, 284488, 402648,
                                                                       178978, 179203, 291208,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 411748, 0, 3,
                                                                       402648, 284698, 402928,
                                                                       179203, 179428, 291523,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 412168, 0, 3,
                                                                       402928, 284908, 403208,
                                                                       179428, 179653, 291838,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 412588, 0, 3,
                                                                       403208, 285118, 403488,
                                                                       179653, 179878, 292153,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 413008, 0, 3,
                                                                       403768, 285538, 404188,
                                                                       180328, 180643, 292468,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 413596, 0, 3,
                                                                       404188, 285853, 404608,
                                                                       180643, 180958, 292909,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 414184, 0, 3,
                                                                       404608, 286168, 405028,
                                                                       180958, 181273, 293350,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 414772, 0, 3,
                                                                       405028, 286483, 405448,
                                                                       181273, 181588, 293791,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 415360, 0, 3,
                                                                       405448, 286798, 405868,
                                                                       181588, 181903, 294232,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 415948, 0, 3,
                                                                       405868, 287113, 406288,
                                                                       181903, 182218, 294673,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 416536, 0, 3,
                                                                       406288, 287428, 406708,
                                                                       182218, 182533, 295114,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 417124, 0, 3,
                                                                       406708, 287743, 407128,
                                                                       182533, 182848, 295555,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 417712, 0, 3,
                                                                       407128, 288058, 407548,
                                                                       182848, 183163, 295996,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 418300, 0, 3,
                                                                       407548, 288373, 407968,
                                                                       183163, 183478, 296437,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 418888, 0, 3,
                                                                       408388, 289003, 408808,
                                                                       184108, 184423, 296878,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 419476, 0, 3,
                                                                       408808, 289318, 409228,
                                                                       184423, 184738, 297319,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 420064, 0, 3,
                                                                       409228, 289633, 409648,
                                                                       184738, 185053, 297760,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 420652, 0, 3,
                                                                       409648, 289948, 410068,
                                                                       185053, 185368, 298201,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 421240, 0, 3,
                                                                       410068, 290263, 410488,
                                                                       185368, 185683, 298642,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 421828, 0, 3,
                                                                       410488, 290578, 410908,
                                                                       185683, 185998, 299083,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 422416, 0, 3,
                                                                       410908, 290893, 411328,
                                                                       185998, 186313, 299524,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 423004, 0, 3,
                                                                       411328, 291208, 411748,
                                                                       186313, 186628, 299965,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 423592, 0, 3,
                                                                       411748, 291523, 412168,
                                                                       186628, 186943, 300406,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 424180, 0, 3,
                                                                       412168, 291838, 412588,
                                                                       186943, 187258, 300847,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 424768, 0, 3,
                                                                       413008, 292468, 413596,
                                                                       187888, 188308, 301288,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 425552, 0, 3,
                                                                       413596, 292909, 414184,
                                                                       188308, 188728, 301876,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 426336, 0, 3,
                                                                       414184, 293350, 414772,
                                                                       188728, 189148, 302464,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 427120, 0, 3,
                                                                       414772, 293791, 415360,
                                                                       189148, 189568, 303052,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 427904, 0, 3,
                                                                       415360, 294232, 415948,
                                                                       189568, 189988, 303640,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 428688, 0, 3,
                                                                       415948, 294673, 416536,
                                                                       189988, 190408, 304228,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 429472, 0, 3,
                                                                       416536, 295114, 417124,
                                                                       190408, 190828, 304816,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 430256, 0, 3,
                                                                       417124, 295555, 417712,
                                                                       190828, 191248, 305404,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 431040, 0, 3,
                                                                       417712, 295996, 418300,
                                                                       191248, 191668, 305992,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 431824, 0, 3,
                                                                       418888, 296878, 419476,
                                                                       192508, 192928, 306580,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 432608, 0, 3,
                                                                       419476, 297319, 420064,
                                                                       192928, 193348, 307168,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 433392, 0, 3,
                                                                       420064, 297760, 420652,
                                                                       193348, 193768, 307756,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 434176, 0, 3,
                                                                       420652, 298201, 421240,
                                                                       193768, 194188, 308344,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 434960, 0, 3,
                                                                       421240, 298642, 421828,
                                                                       194188, 194608, 308932,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 435744, 0, 3,
                                                                       421828, 299083, 422416,
                                                                       194608, 195028, 309520,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 436528, 0, 3,
                                                                       422416, 299524, 423004,
                                                                       195028, 195448, 310108,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 437312, 0, 3,
                                                                       423004, 299965, 423592,
                                                                       195448, 195868, 310696,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 438096, 0, 3,
                                                                       423592, 300406, 424180,
                                                                       195868, 196288, 311284,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 438880, 0, 3,
                                                                       424768, 301288, 425552,
                                                                       197128, 197668, 311872,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 439888, 0, 3,
                                                                       425552, 301876, 426336,
                                                                       197668, 198208, 312628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 440896, 0, 3,
                                                                       426336, 302464, 427120,
                                                                       198208, 198748, 313384,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 441904, 0, 3,
                                                                       427120, 303052, 427904,
                                                                       198748, 199288, 314140,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 442912, 0, 3,
                                                                       427904, 303640, 428688,
                                                                       199288, 199828, 314896,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 443920, 0, 3,
                                                                       428688, 304228, 429472,
                                                                       199828, 200368, 315652,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 444928, 0, 3,
                                                                       429472, 304816, 430256,
                                                                       200368, 200908, 316408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 445936, 0, 3,
                                                                       430256, 305404, 431040,
                                                                       200908, 201448, 317164,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 446944, 0, 3,
                                                                       431824, 306580, 432608,
                                                                       202528, 203068, 317920,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 447952, 0, 3,
                                                                       432608, 307168, 433392,
                                                                       203068, 203608, 318676,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 448960, 0, 3,
                                                                       433392, 307756, 434176,
                                                                       203608, 204148, 319432,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 449968, 0, 3,
                                                                       434176, 308344, 434960,
                                                                       204148, 204688, 320188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 450976, 0, 3,
                                                                       434960, 308932, 435744,
                                                                       204688, 205228, 320944,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 451984, 0, 3,
                                                                       435744, 309520, 436528,
                                                                       205228, 205768, 321700,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 452992, 0, 3,
                                                                       436528, 310108, 437312,
                                                                       205768, 206308, 322456,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 454000, 0, 3,
                                                                       437312, 310696, 438096,
                                                                       206308, 206848, 323212,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 455008, 0, 3,
                                                                       438880, 311872, 439888,
                                                                       207928, 208603, 323968,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 456268, 0, 3,
                                                                       439888, 312628, 440896,
                                                                       208603, 209278, 324913,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 457528, 0, 3,
                                                                       440896, 313384, 441904,
                                                                       209278, 209953, 325858,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 458788, 0, 3,
                                                                       441904, 314140, 442912,
                                                                       209953, 210628, 326803,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 460048, 0, 3,
                                                                       442912, 314896, 443920,
                                                                       210628, 211303, 327748,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 461308, 0, 3,
                                                                       443920, 315652, 444928,
                                                                       211303, 211978, 328693,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 462568, 0, 3,
                                                                       444928, 316408, 445936,
                                                                       211978, 212653, 329638,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 463828, 0, 3,
                                                                       446944, 317920, 447952,
                                                                       214003, 214678, 330583,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 465088, 0, 3,
                                                                       447952, 318676, 448960,
                                                                       214678, 215353, 331528,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 466348, 0, 3,
                                                                       448960, 319432, 449968,
                                                                       215353, 216028, 332473,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 467608, 0, 3,
                                                                       449968, 320188, 450976,
                                                                       216028, 216703, 333418,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 468868, 0, 3,
                                                                       450976, 320944, 451984,
                                                                       216703, 217378, 334363,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 470128, 0, 3,
                                                                       451984, 321700, 452992,
                                                                       217378, 218053, 335308,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 471388, 0, 3,
                                                                       452992, 322456, 454000,
                                                                       218053, 218728, 336253,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 472648, 0, 3,
                                                                       455008, 323968, 456268,
                                                                       220078, 220903, 337198,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 474188, 0, 3,
                                                                       456268, 324913, 457528,
                                                                       220903, 221728, 338353,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 475728, 0, 3,
                                                                       457528, 325858, 458788,
                                                                       221728, 222553, 339508,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 477268, 0, 3,
                                                                       458788, 326803, 460048,
                                                                       222553, 223378, 340663,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 478808, 0, 3,
                                                                       460048, 327748, 461308,
                                                                       223378, 224203, 341818,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 480348, 0, 3,
                                                                       461308, 328693, 462568,
                                                                       224203, 225028, 342973,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 481888, 0, 3,
                                                                       463828, 330583, 465088,
                                                                       226678, 227503, 344128,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 483428, 0, 3,
                                                                       465088, 331528, 466348,
                                                                       227503, 228328, 345283,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 484968, 0, 3,
                                                                       466348, 332473, 467608,
                                                                       228328, 229153, 346438,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 486508, 0, 3,
                                                                       467608, 333418, 468868,
                                                                       229153, 229978, 347593,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 488048, 0, 3,
                                                                       468868, 334363, 470128,
                                                                       229978, 230803, 348748,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 489588, 0, 3,
                                                                       470128, 335308, 471388,
                                                                       230803, 231628, 349903,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 491128, 0, 3,
                                                                       472648, 337198, 474188,
                                                                       233278, 234268, 351058,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 492976, 0, 3,
                                                                       474188, 338353, 475728,
                                                                       234268, 235258, 352444,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 494824, 0, 3,
                                                                       475728, 339508, 477268,
                                                                       235258, 236248, 353830,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 496672, 0, 3,
                                                                       477268, 340663, 478808,
                                                                       236248, 237238, 355216,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 498520, 0, 3,
                                                                       478808, 341818, 480348,
                                                                       237238, 238228, 356602,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 500368, 0, 3,
                                                                       481888, 344128, 483428,
                                                                       240208, 241198, 357988,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 502216, 0, 3,
                                                                       483428, 345283, 484968,
                                                                       241198, 242188, 359374,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 504064, 0, 3,
                                                                       484968, 346438, 486508,
                                                                       242188, 243178, 360760,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 505912, 0, 3,
                                                                       486508, 347593, 488048,
                                                                       243178, 244168, 362146,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 507760, 0, 3,
                                                                       488048, 348748, 489588,
                                                                       244168, 245158, 363532,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 509608, 0, 3,
                                                                       491128, 351058, 492976,
                                                                       247138, 248308, 364918,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 511792, 0, 3,
                                                                       492976, 352444, 494824,
                                                                       248308, 249478, 366556,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 513976, 0, 3,
                                                                       494824, 353830, 496672,
                                                                       249478, 250648, 368194,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 516160, 0, 3,
                                                                       496672, 355216, 498520,
                                                                       250648, 251818, 369832,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 518344, 0, 3,
                                                                       500368, 357988, 502216,
                                                                       254158, 255328, 371470,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 520528, 0, 3,
                                                                       502216, 359374, 504064,
                                                                       255328, 256498, 373108,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 522712, 0, 3,
                                                                       504064, 360760, 505912,
                                                                       256498, 257668, 374746,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 524896, 0, 3,
                                                                       505912, 362146, 507760,
                                                                       257668, 258838, 376384,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 527080, 0, 3,
                                                                       509608, 364918, 511792,
                                                                       261178, 262543, 378022,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 529628, 0, 3,
                                                                       511792, 366556, 513976,
                                                                       262543, 263908, 379933,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 532176, 0, 3,
                                                                       513976, 368194, 516160,
                                                                       263908, 265273, 381844,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 534724, 0, 3,
                                                                       518344, 371470, 520528,
                                                                       268003, 269368, 383755,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 537272, 0, 3,
                                                                       520528, 373108, 522712,
                                                                       269368, 270733, 385666,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 539820, 0, 3,
                                                                       522712, 374746, 524896,
                                                                       270733, 272098, 387577,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542368, 3, 274828,
                                                                       274849, 389544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542404, 3, 274849,
                                                                       274870, 389572, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542440, 3, 274870,
                                                                       274891, 389600, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542476, 3, 274891,
                                                                       274912, 389628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542512, 3, 274912,
                                                                       274933, 389656, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542548, 3, 274933,
                                                                       274954, 389684, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542584, 3, 274954,
                                                                       274975, 389712, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542620, 3, 274975,
                                                                       274996, 389740, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542656, 3, 274996,
                                                                       275017, 389768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542692, 3, 275017,
                                                                       275038, 389796, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542728, 3, 275038,
                                                                       275059, 389824, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542764, 3, 275059,
                                                                       275080, 389852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542800, 3, 275080,
                                                                       275101, 389880, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542836, 3, 275143,
                                                                       275164, 389964, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542872, 3, 275164,
                                                                       275185, 389992, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542908, 3, 275185,
                                                                       275206, 390020, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542944, 3, 275206,
                                                                       275227, 390048, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 542980, 3, 275227,
                                                                       275248, 390076, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 543016, 3, 275248,
                                                                       275269, 390104, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 543052, 3, 275269,
                                                                       275290, 390132, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 543088, 3, 275290,
                                                                       275311, 390160, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 543124, 3, 275311,
                                                                       275332, 390188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 543160, 3, 275332,
                                                                       275353, 390216, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 543196, 3, 275353,
                                                                       275374, 390244, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 543232, 3, 275374,
                                                                       275395, 390272, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 543268, 3, 275395,
                                                                       275416, 390300, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 543304, 0, 3,
                                                                       542368, 389544, 542404,
                                                                       275458, 275521, 390496,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 543412, 0, 3,
                                                                       542404, 389572, 542440,
                                                                       275521, 275584, 390580,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 543520, 0, 3,
                                                                       542440, 389600, 542476,
                                                                       275584, 275647, 390664,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 543628, 0, 3,
                                                                       542476, 389628, 542512,
                                                                       275647, 275710, 390748,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 543736, 0, 3,
                                                                       542512, 389656, 542548,
                                                                       275710, 275773, 390832,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 543844, 0, 3,
                                                                       542548, 389684, 542584,
                                                                       275773, 275836, 390916,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 543952, 0, 3,
                                                                       542584, 389712, 542620,
                                                                       275836, 275899, 391000,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 544060, 0, 3,
                                                                       542620, 389740, 542656,
                                                                       275899, 275962, 391084,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 544168, 0, 3,
                                                                       542656, 389768, 542692,
                                                                       275962, 276025, 391168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 544276, 0, 3,
                                                                       542692, 389796, 542728,
                                                                       276025, 276088, 391252,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 544384, 0, 3,
                                                                       542728, 389824, 542764,
                                                                       276088, 276151, 391336,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 544492, 0, 3,
                                                                       542764, 389852, 542800,
                                                                       276151, 276214, 391420,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 544600, 0, 3,
                                                                       542836, 389964, 542872,
                                                                       276340, 276403, 391672,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 544708, 0, 3,
                                                                       542872, 389992, 542908,
                                                                       276403, 276466, 391756,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 544816, 0, 3,
                                                                       542908, 390020, 542944,
                                                                       276466, 276529, 391840,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 544924, 0, 3,
                                                                       542944, 390048, 542980,
                                                                       276529, 276592, 391924,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 545032, 0, 3,
                                                                       542980, 390076, 543016,
                                                                       276592, 276655, 392008,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 545140, 0, 3,
                                                                       543016, 390104, 543052,
                                                                       276655, 276718, 392092,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 545248, 0, 3,
                                                                       543052, 390132, 543088,
                                                                       276718, 276781, 392176,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 545356, 0, 3,
                                                                       543088, 390160, 543124,
                                                                       276781, 276844, 392260,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 545464, 0, 3,
                                                                       543124, 390188, 543160,
                                                                       276844, 276907, 392344,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 545572, 0, 3,
                                                                       543160, 390216, 543196,
                                                                       276907, 276970, 392428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 545680, 0, 3,
                                                                       543196, 390244, 543232,
                                                                       276970, 277033, 392512,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 545788, 0, 3,
                                                                       543232, 390272, 543268,
                                                                       277033, 277096, 392596,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 545896, 0, 3,
                                                                       543304, 390496, 543412,
                                                                       277222, 277348, 393016,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 546112, 0, 3,
                                                                       543412, 390580, 543520,
                                                                       277348, 277474, 393184,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 546328, 0, 3,
                                                                       543520, 390664, 543628,
                                                                       277474, 277600, 393352,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 546544, 0, 3,
                                                                       543628, 390748, 543736,
                                                                       277600, 277726, 393520,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 546760, 0, 3,
                                                                       543736, 390832, 543844,
                                                                       277726, 277852, 393688,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 546976, 0, 3,
                                                                       543844, 390916, 543952,
                                                                       277852, 277978, 393856,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 547192, 0, 3,
                                                                       543952, 391000, 544060,
                                                                       277978, 278104, 394024,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 547408, 0, 3,
                                                                       544060, 391084, 544168,
                                                                       278104, 278230, 394192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 547624, 0, 3,
                                                                       544168, 391168, 544276,
                                                                       278230, 278356, 394360,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 547840, 0, 3,
                                                                       544276, 391252, 544384,
                                                                       278356, 278482, 394528,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 548056, 0, 3,
                                                                       544384, 391336, 544492,
                                                                       278482, 278608, 394696,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 548272, 0, 3,
                                                                       544600, 391672, 544708,
                                                                       278860, 278986, 395200,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 548488, 0, 3,
                                                                       544708, 391756, 544816,
                                                                       278986, 279112, 395368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 548704, 0, 3,
                                                                       544816, 391840, 544924,
                                                                       279112, 279238, 395536,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 548920, 0, 3,
                                                                       544924, 391924, 545032,
                                                                       279238, 279364, 395704,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 549136, 0, 3,
                                                                       545032, 392008, 545140,
                                                                       279364, 279490, 395872,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 549352, 0, 3,
                                                                       545140, 392092, 545248,
                                                                       279490, 279616, 396040,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 549568, 0, 3,
                                                                       545248, 392176, 545356,
                                                                       279616, 279742, 396208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 549784, 0, 3,
                                                                       545356, 392260, 545464,
                                                                       279742, 279868, 396376,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 550000, 0, 3,
                                                                       545464, 392344, 545572,
                                                                       279868, 279994, 396544,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 550216, 0, 3,
                                                                       545572, 392428, 545680,
                                                                       279994, 280120, 396712,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 550432, 0, 3,
                                                                       545680, 392512, 545788,
                                                                       280120, 280246, 396880,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 550648, 0, 3,
                                                                       545896, 393016, 546112,
                                                                       280498, 280708, 397608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 551008, 0, 3,
                                                                       546112, 393184, 546328,
                                                                       280708, 280918, 397888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 551368, 0, 3,
                                                                       546328, 393352, 546544,
                                                                       280918, 281128, 398168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 551728, 0, 3,
                                                                       546544, 393520, 546760,
                                                                       281128, 281338, 398448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 552088, 0, 3,
                                                                       546760, 393688, 546976,
                                                                       281338, 281548, 398728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 552448, 0, 3,
                                                                       546976, 393856, 547192,
                                                                       281548, 281758, 399008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 552808, 0, 3,
                                                                       547192, 394024, 547408,
                                                                       281758, 281968, 399288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 553168, 0, 3,
                                                                       547408, 394192, 547624,
                                                                       281968, 282178, 399568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 553528, 0, 3,
                                                                       547624, 394360, 547840,
                                                                       282178, 282388, 399848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 553888, 0, 3,
                                                                       547840, 394528, 548056,
                                                                       282388, 282598, 400128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 554248, 0, 3,
                                                                       548272, 395200, 548488,
                                                                       283018, 283228, 400968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 554608, 0, 3,
                                                                       548488, 395368, 548704,
                                                                       283228, 283438, 401248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 554968, 0, 3,
                                                                       548704, 395536, 548920,
                                                                       283438, 283648, 401528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 555328, 0, 3,
                                                                       548920, 395704, 549136,
                                                                       283648, 283858, 401808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 555688, 0, 3,
                                                                       549136, 395872, 549352,
                                                                       283858, 284068, 402088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 556048, 0, 3,
                                                                       549352, 396040, 549568,
                                                                       284068, 284278, 402368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 556408, 0, 3,
                                                                       549568, 396208, 549784,
                                                                       284278, 284488, 402648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 556768, 0, 3,
                                                                       549784, 396376, 550000,
                                                                       284488, 284698, 402928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 557128, 0, 3,
                                                                       550000, 396544, 550216,
                                                                       284698, 284908, 403208,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 557488, 0, 3,
                                                                       550216, 396712, 550432,
                                                                       284908, 285118, 403488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 557848, 0, 3,
                                                                       550648, 397608, 551008,
                                                                       285538, 285853, 404608,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 558388, 0, 3,
                                                                       551008, 397888, 551368,
                                                                       285853, 286168, 405028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 558928, 0, 3,
                                                                       551368, 398168, 551728,
                                                                       286168, 286483, 405448,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 559468, 0, 3,
                                                                       551728, 398448, 552088,
                                                                       286483, 286798, 405868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 560008, 0, 3,
                                                                       552088, 398728, 552448,
                                                                       286798, 287113, 406288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 560548, 0, 3,
                                                                       552448, 399008, 552808,
                                                                       287113, 287428, 406708,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 561088, 0, 3,
                                                                       552808, 399288, 553168,
                                                                       287428, 287743, 407128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 561628, 0, 3,
                                                                       553168, 399568, 553528,
                                                                       287743, 288058, 407548,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 562168, 0, 3,
                                                                       553528, 399848, 553888,
                                                                       288058, 288373, 407968,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 562708, 0, 3,
                                                                       554248, 400968, 554608,
                                                                       289003, 289318, 409228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 563248, 0, 3,
                                                                       554608, 401248, 554968,
                                                                       289318, 289633, 409648,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 563788, 0, 3,
                                                                       554968, 401528, 555328,
                                                                       289633, 289948, 410068,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 564328, 0, 3,
                                                                       555328, 401808, 555688,
                                                                       289948, 290263, 410488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 564868, 0, 3,
                                                                       555688, 402088, 556048,
                                                                       290263, 290578, 410908,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 565408, 0, 3,
                                                                       556048, 402368, 556408,
                                                                       290578, 290893, 411328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 565948, 0, 3,
                                                                       556408, 402648, 556768,
                                                                       290893, 291208, 411748,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 566488, 0, 3,
                                                                       556768, 402928, 557128,
                                                                       291208, 291523, 412168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 567028, 0, 3,
                                                                       557128, 403208, 557488,
                                                                       291523, 291838, 412588,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 567568, 0, 3,
                                                                       557848, 404608, 558388,
                                                                       292468, 292909, 414184,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 568324, 0, 3,
                                                                       558388, 405028, 558928,
                                                                       292909, 293350, 414772,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 569080, 0, 3,
                                                                       558928, 405448, 559468,
                                                                       293350, 293791, 415360,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 569836, 0, 3,
                                                                       559468, 405868, 560008,
                                                                       293791, 294232, 415948,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 570592, 0, 3,
                                                                       560008, 406288, 560548,
                                                                       294232, 294673, 416536,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 571348, 0, 3,
                                                                       560548, 406708, 561088,
                                                                       294673, 295114, 417124,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 572104, 0, 3,
                                                                       561088, 407128, 561628,
                                                                       295114, 295555, 417712,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 572860, 0, 3,
                                                                       561628, 407548, 562168,
                                                                       295555, 295996, 418300,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 573616, 0, 3,
                                                                       562708, 409228, 563248,
                                                                       296878, 297319, 420064,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 574372, 0, 3,
                                                                       563248, 409648, 563788,
                                                                       297319, 297760, 420652,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 575128, 0, 3,
                                                                       563788, 410068, 564328,
                                                                       297760, 298201, 421240,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 575884, 0, 3,
                                                                       564328, 410488, 564868,
                                                                       298201, 298642, 421828,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 576640, 0, 3,
                                                                       564868, 410908, 565408,
                                                                       298642, 299083, 422416,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 577396, 0, 3,
                                                                       565408, 411328, 565948,
                                                                       299083, 299524, 423004,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 578152, 0, 3,
                                                                       565948, 411748, 566488,
                                                                       299524, 299965, 423592,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 578908, 0, 3,
                                                                       566488, 412168, 567028,
                                                                       299965, 300406, 424180,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 579664, 0, 3,
                                                                       567568, 414184, 568324,
                                                                       301288, 301876, 426336,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 580672, 0, 3,
                                                                       568324, 414772, 569080,
                                                                       301876, 302464, 427120,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 581680, 0, 3,
                                                                       569080, 415360, 569836,
                                                                       302464, 303052, 427904,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 582688, 0, 3,
                                                                       569836, 415948, 570592,
                                                                       303052, 303640, 428688,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 583696, 0, 3,
                                                                       570592, 416536, 571348,
                                                                       303640, 304228, 429472,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 584704, 0, 3,
                                                                       571348, 417124, 572104,
                                                                       304228, 304816, 430256,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 585712, 0, 3,
                                                                       572104, 417712, 572860,
                                                                       304816, 305404, 431040,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 586720, 0, 3,
                                                                       573616, 420064, 574372,
                                                                       306580, 307168, 433392,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 587728, 0, 3,
                                                                       574372, 420652, 575128,
                                                                       307168, 307756, 434176,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 588736, 0, 3,
                                                                       575128, 421240, 575884,
                                                                       307756, 308344, 434960,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 589744, 0, 3,
                                                                       575884, 421828, 576640,
                                                                       308344, 308932, 435744,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 590752, 0, 3,
                                                                       576640, 422416, 577396,
                                                                       308932, 309520, 436528,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 591760, 0, 3,
                                                                       577396, 423004, 578152,
                                                                       309520, 310108, 437312,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 592768, 0, 3,
                                                                       578152, 423592, 578908,
                                                                       310108, 310696, 438096,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 593776, 0, 3,
                                                                       579664, 426336, 580672,
                                                                       311872, 312628, 440896,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 595072, 0, 3,
                                                                       580672, 427120, 581680,
                                                                       312628, 313384, 441904,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 596368, 0, 3,
                                                                       581680, 427904, 582688,
                                                                       313384, 314140, 442912,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 597664, 0, 3,
                                                                       582688, 428688, 583696,
                                                                       314140, 314896, 443920,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 598960, 0, 3,
                                                                       583696, 429472, 584704,
                                                                       314896, 315652, 444928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 600256, 0, 3,
                                                                       584704, 430256, 585712,
                                                                       315652, 316408, 445936,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 601552, 0, 3,
                                                                       586720, 433392, 587728,
                                                                       317920, 318676, 448960,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 602848, 0, 3,
                                                                       587728, 434176, 588736,
                                                                       318676, 319432, 449968,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 604144, 0, 3,
                                                                       588736, 434960, 589744,
                                                                       319432, 320188, 450976,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 605440, 0, 3,
                                                                       589744, 435744, 590752,
                                                                       320188, 320944, 451984,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 606736, 0, 3,
                                                                       590752, 436528, 591760,
                                                                       320944, 321700, 452992,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 608032, 0, 3,
                                                                       591760, 437312, 592768,
                                                                       321700, 322456, 454000,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 609328, 0, 3,
                                                                       593776, 440896, 595072,
                                                                       323968, 324913, 457528,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 610948, 0, 3,
                                                                       595072, 441904, 596368,
                                                                       324913, 325858, 458788,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 612568, 0, 3,
                                                                       596368, 442912, 597664,
                                                                       325858, 326803, 460048,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 614188, 0, 3,
                                                                       597664, 443920, 598960,
                                                                       326803, 327748, 461308,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 615808, 0, 3,
                                                                       598960, 444928, 600256,
                                                                       327748, 328693, 462568,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 617428, 0, 3,
                                                                       601552, 448960, 602848,
                                                                       330583, 331528, 466348,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 619048, 0, 3,
                                                                       602848, 449968, 604144,
                                                                       331528, 332473, 467608,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 620668, 0, 3,
                                                                       604144, 450976, 605440,
                                                                       332473, 333418, 468868,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 622288, 0, 3,
                                                                       605440, 451984, 606736,
                                                                       333418, 334363, 470128,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 623908, 0, 3,
                                                                       606736, 452992, 608032,
                                                                       334363, 335308, 471388,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 625528, 0, 3,
                                                                       609328, 457528, 610948,
                                                                       337198, 338353, 475728,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 627508, 0, 3,
                                                                       610948, 458788, 612568,
                                                                       338353, 339508, 477268,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 629488, 0, 3,
                                                                       612568, 460048, 614188,
                                                                       339508, 340663, 478808,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 631468, 0, 3,
                                                                       614188, 461308, 615808,
                                                                       340663, 341818, 480348,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 633448, 0, 3,
                                                                       617428, 466348, 619048,
                                                                       344128, 345283, 484968,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 635428, 0, 3,
                                                                       619048, 467608, 620668,
                                                                       345283, 346438, 486508,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 637408, 0, 3,
                                                                       620668, 468868, 622288,
                                                                       346438, 347593, 488048,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 639388, 0, 3,
                                                                       622288, 470128, 623908,
                                                                       347593, 348748, 489588,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 641368, 0, 3,
                                                                       625528, 475728, 627508,
                                                                       351058, 352444, 494824,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 643744, 0, 3,
                                                                       627508, 477268, 629488,
                                                                       352444, 353830, 496672,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 646120, 0, 3,
                                                                       629488, 478808, 631468,
                                                                       353830, 355216, 498520,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 648496, 0, 3,
                                                                       633448, 484968, 635428,
                                                                       357988, 359374, 504064,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 650872, 0, 3,
                                                                       635428, 486508, 637408,
                                                                       359374, 360760, 505912,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 653248, 0, 3,
                                                                       637408, 488048, 639388,
                                                                       360760, 362146, 507760,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 655624, 0, 3,
                                                                       641368, 494824, 643744,
                                                                       364918, 366556, 513976,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 658432, 0, 3,
                                                                       643744, 496672, 646120,
                                                                       366556, 368194, 516160,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 661240, 0, 3,
                                                                       648496, 504064, 650872,
                                                                       371470, 373108, 522712,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 664048, 0, 3,
                                                                       650872, 505912, 653248,
                                                                       373108, 374746, 524896,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsk_three_center_electron_repulsion_0(buffer, 666856, 0, 3,
                                                                       655624, 513976, 658432,
                                                                       378022, 379933, 532176,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsk_three_center_electron_repulsion_0(buffer, 670132, 0, 3,
                                                                       661240, 522712, 664048,
                                                                       383755, 385666, 539820,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673408, 3, 389488,
                                                                       389516, 542368, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673453, 3, 389516,
                                                                       389544, 542404, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673498, 3, 389544,
                                                                       389572, 542440, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673543, 3, 389572,
                                                                       389600, 542476, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673588, 3, 389600,
                                                                       389628, 542512, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673633, 3, 389628,
                                                                       389656, 542548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673678, 3, 389656,
                                                                       389684, 542584, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673723, 3, 389684,
                                                                       389712, 542620, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673768, 3, 389712,
                                                                       389740, 542656, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673813, 3, 389740,
                                                                       389768, 542692, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673858, 3, 389768,
                                                                       389796, 542728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673903, 3, 389796,
                                                                       389824, 542764, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673948, 3, 389824,
                                                                       389852, 542800, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 673993, 3, 389908,
                                                                       389936, 542836, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674038, 3, 389936,
                                                                       389964, 542872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674083, 3, 389964,
                                                                       389992, 542908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674128, 3, 389992,
                                                                       390020, 542944, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674173, 3, 390020,
                                                                       390048, 542980, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674218, 3, 390048,
                                                                       390076, 543016, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674263, 3, 390076,
                                                                       390104, 543052, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674308, 3, 390104,
                                                                       390132, 543088, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674353, 3, 390132,
                                                                       390160, 543124, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674398, 3, 390160,
                                                                       390188, 543160, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674443, 3, 390188,
                                                                       390216, 543196, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674488, 3, 390216,
                                                                       390244, 543232, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 674533, 3, 390244,
                                                                       390272, 543268, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 674578, 0, 3,
                                                                       673408, 542368, 673453,
                                                                       390328, 390412, 543304,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 674713, 0, 3,
                                                                       673453, 542404, 673498,
                                                                       390412, 390496, 543412,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 674848, 0, 3,
                                                                       673498, 542440, 673543,
                                                                       390496, 390580, 543520,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 674983, 0, 3,
                                                                       673543, 542476, 673588,
                                                                       390580, 390664, 543628,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 675118, 0, 3,
                                                                       673588, 542512, 673633,
                                                                       390664, 390748, 543736,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 675253, 0, 3,
                                                                       673633, 542548, 673678,
                                                                       390748, 390832, 543844,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 675388, 0, 3,
                                                                       673678, 542584, 673723,
                                                                       390832, 390916, 543952,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 675523, 0, 3,
                                                                       673723, 542620, 673768,
                                                                       390916, 391000, 544060,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 675658, 0, 3,
                                                                       673768, 542656, 673813,
                                                                       391000, 391084, 544168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 675793, 0, 3,
                                                                       673813, 542692, 673858,
                                                                       391084, 391168, 544276,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 675928, 0, 3,
                                                                       673858, 542728, 673903,
                                                                       391168, 391252, 544384,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 676063, 0, 3,
                                                                       673903, 542764, 673948,
                                                                       391252, 391336, 544492,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 676198, 0, 3,
                                                                       673993, 542836, 674038,
                                                                       391504, 391588, 544600,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 676333, 0, 3,
                                                                       674038, 542872, 674083,
                                                                       391588, 391672, 544708,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 676468, 0, 3,
                                                                       674083, 542908, 674128,
                                                                       391672, 391756, 544816,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 676603, 0, 3,
                                                                       674128, 542944, 674173,
                                                                       391756, 391840, 544924,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 676738, 0, 3,
                                                                       674173, 542980, 674218,
                                                                       391840, 391924, 545032,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 676873, 0, 3,
                                                                       674218, 543016, 674263,
                                                                       391924, 392008, 545140,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 677008, 0, 3,
                                                                       674263, 543052, 674308,
                                                                       392008, 392092, 545248,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 677143, 0, 3,
                                                                       674308, 543088, 674353,
                                                                       392092, 392176, 545356,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 677278, 0, 3,
                                                                       674353, 543124, 674398,
                                                                       392176, 392260, 545464,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 677413, 0, 3,
                                                                       674398, 543160, 674443,
                                                                       392260, 392344, 545572,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 677548, 0, 3,
                                                                       674443, 543196, 674488,
                                                                       392344, 392428, 545680,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 677683, 0, 3,
                                                                       674488, 543232, 674533,
                                                                       392428, 392512, 545788,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 677818, 0, 3,
                                                                       674578, 543304, 674713,
                                                                       392680, 392848, 545896,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 678088, 0, 3,
                                                                       674713, 543412, 674848,
                                                                       392848, 393016, 546112,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 678358, 0, 3,
                                                                       674848, 543520, 674983,
                                                                       393016, 393184, 546328,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 678628, 0, 3,
                                                                       674983, 543628, 675118,
                                                                       393184, 393352, 546544,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 678898, 0, 3,
                                                                       675118, 543736, 675253,
                                                                       393352, 393520, 546760,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 679168, 0, 3,
                                                                       675253, 543844, 675388,
                                                                       393520, 393688, 546976,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 679438, 0, 3,
                                                                       675388, 543952, 675523,
                                                                       393688, 393856, 547192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 679708, 0, 3,
                                                                       675523, 544060, 675658,
                                                                       393856, 394024, 547408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 679978, 0, 3,
                                                                       675658, 544168, 675793,
                                                                       394024, 394192, 547624,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 680248, 0, 3,
                                                                       675793, 544276, 675928,
                                                                       394192, 394360, 547840,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 680518, 0, 3,
                                                                       675928, 544384, 676063,
                                                                       394360, 394528, 548056,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 680788, 0, 3,
                                                                       676198, 544600, 676333,
                                                                       394864, 395032, 548272,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 681058, 0, 3,
                                                                       676333, 544708, 676468,
                                                                       395032, 395200, 548488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 681328, 0, 3,
                                                                       676468, 544816, 676603,
                                                                       395200, 395368, 548704,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 681598, 0, 3,
                                                                       676603, 544924, 676738,
                                                                       395368, 395536, 548920,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 681868, 0, 3,
                                                                       676738, 545032, 676873,
                                                                       395536, 395704, 549136,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 682138, 0, 3,
                                                                       676873, 545140, 677008,
                                                                       395704, 395872, 549352,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 682408, 0, 3,
                                                                       677008, 545248, 677143,
                                                                       395872, 396040, 549568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 682678, 0, 3,
                                                                       677143, 545356, 677278,
                                                                       396040, 396208, 549784,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 682948, 0, 3,
                                                                       677278, 545464, 677413,
                                                                       396208, 396376, 550000,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 683218, 0, 3,
                                                                       677413, 545572, 677548,
                                                                       396376, 396544, 550216,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 683488, 0, 3,
                                                                       677548, 545680, 677683,
                                                                       396544, 396712, 550432,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 683758, 0, 3,
                                                                       677818, 545896, 678088,
                                                                       397048, 397328, 550648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 684208, 0, 3,
                                                                       678088, 546112, 678358,
                                                                       397328, 397608, 551008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 684658, 0, 3,
                                                                       678358, 546328, 678628,
                                                                       397608, 397888, 551368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 685108, 0, 3,
                                                                       678628, 546544, 678898,
                                                                       397888, 398168, 551728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 685558, 0, 3,
                                                                       678898, 546760, 679168,
                                                                       398168, 398448, 552088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 686008, 0, 3,
                                                                       679168, 546976, 679438,
                                                                       398448, 398728, 552448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 686458, 0, 3,
                                                                       679438, 547192, 679708,
                                                                       398728, 399008, 552808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 686908, 0, 3,
                                                                       679708, 547408, 679978,
                                                                       399008, 399288, 553168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 687358, 0, 3,
                                                                       679978, 547624, 680248,
                                                                       399288, 399568, 553528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 687808, 0, 3,
                                                                       680248, 547840, 680518,
                                                                       399568, 399848, 553888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 688258, 0, 3,
                                                                       680788, 548272, 681058,
                                                                       400408, 400688, 554248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 688708, 0, 3,
                                                                       681058, 548488, 681328,
                                                                       400688, 400968, 554608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 689158, 0, 3,
                                                                       681328, 548704, 681598,
                                                                       400968, 401248, 554968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 689608, 0, 3,
                                                                       681598, 548920, 681868,
                                                                       401248, 401528, 555328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 690058, 0, 3,
                                                                       681868, 549136, 682138,
                                                                       401528, 401808, 555688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 690508, 0, 3,
                                                                       682138, 549352, 682408,
                                                                       401808, 402088, 556048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 690958, 0, 3,
                                                                       682408, 549568, 682678,
                                                                       402088, 402368, 556408,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 691408, 0, 3,
                                                                       682678, 549784, 682948,
                                                                       402368, 402648, 556768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 691858, 0, 3,
                                                                       682948, 550000, 683218,
                                                                       402648, 402928, 557128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 692308, 0, 3,
                                                                       683218, 550216, 683488,
                                                                       402928, 403208, 557488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 692758, 0, 3,
                                                                       683758, 550648, 684208,
                                                                       403768, 404188, 557848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 693433, 0, 3,
                                                                       684208, 551008, 684658,
                                                                       404188, 404608, 558388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 694108, 0, 3,
                                                                       684658, 551368, 685108,
                                                                       404608, 405028, 558928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 694783, 0, 3,
                                                                       685108, 551728, 685558,
                                                                       405028, 405448, 559468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 695458, 0, 3,
                                                                       685558, 552088, 686008,
                                                                       405448, 405868, 560008,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 696133, 0, 3,
                                                                       686008, 552448, 686458,
                                                                       405868, 406288, 560548,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 696808, 0, 3,
                                                                       686458, 552808, 686908,
                                                                       406288, 406708, 561088,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 697483, 0, 3,
                                                                       686908, 553168, 687358,
                                                                       406708, 407128, 561628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 698158, 0, 3,
                                                                       687358, 553528, 687808,
                                                                       407128, 407548, 562168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 698833, 0, 3,
                                                                       688258, 554248, 688708,
                                                                       408388, 408808, 562708,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 699508, 0, 3,
                                                                       688708, 554608, 689158,
                                                                       408808, 409228, 563248,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 700183, 0, 3,
                                                                       689158, 554968, 689608,
                                                                       409228, 409648, 563788,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 700858, 0, 3,
                                                                       689608, 555328, 690058,
                                                                       409648, 410068, 564328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 701533, 0, 3,
                                                                       690058, 555688, 690508,
                                                                       410068, 410488, 564868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 702208, 0, 3,
                                                                       690508, 556048, 690958,
                                                                       410488, 410908, 565408,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 702883, 0, 3,
                                                                       690958, 556408, 691408,
                                                                       410908, 411328, 565948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 703558, 0, 3,
                                                                       691408, 556768, 691858,
                                                                       411328, 411748, 566488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 704233, 0, 3,
                                                                       691858, 557128, 692308,
                                                                       411748, 412168, 567028,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 704908, 0, 3,
                                                                       692758, 557848, 693433,
                                                                       413008, 413596, 567568,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 705853, 0, 3,
                                                                       693433, 558388, 694108,
                                                                       413596, 414184, 568324,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 706798, 0, 3,
                                                                       694108, 558928, 694783,
                                                                       414184, 414772, 569080,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 707743, 0, 3,
                                                                       694783, 559468, 695458,
                                                                       414772, 415360, 569836,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 708688, 0, 3,
                                                                       695458, 560008, 696133,
                                                                       415360, 415948, 570592,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 709633, 0, 3,
                                                                       696133, 560548, 696808,
                                                                       415948, 416536, 571348,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 710578, 0, 3,
                                                                       696808, 561088, 697483,
                                                                       416536, 417124, 572104,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 711523, 0, 3,
                                                                       697483, 561628, 698158,
                                                                       417124, 417712, 572860,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 712468, 0, 3,
                                                                       698833, 562708, 699508,
                                                                       418888, 419476, 573616,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 713413, 0, 3,
                                                                       699508, 563248, 700183,
                                                                       419476, 420064, 574372,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 714358, 0, 3,
                                                                       700183, 563788, 700858,
                                                                       420064, 420652, 575128,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 715303, 0, 3,
                                                                       700858, 564328, 701533,
                                                                       420652, 421240, 575884,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 716248, 0, 3,
                                                                       701533, 564868, 702208,
                                                                       421240, 421828, 576640,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 717193, 0, 3,
                                                                       702208, 565408, 702883,
                                                                       421828, 422416, 577396,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 718138, 0, 3,
                                                                       702883, 565948, 703558,
                                                                       422416, 423004, 578152,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 719083, 0, 3,
                                                                       703558, 566488, 704233,
                                                                       423004, 423592, 578908,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 720028, 0, 3,
                                                                       704908, 567568, 705853,
                                                                       424768, 425552, 579664,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 721288, 0, 3,
                                                                       705853, 568324, 706798,
                                                                       425552, 426336, 580672,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 722548, 0, 3,
                                                                       706798, 569080, 707743,
                                                                       426336, 427120, 581680,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 723808, 0, 3,
                                                                       707743, 569836, 708688,
                                                                       427120, 427904, 582688,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 725068, 0, 3,
                                                                       708688, 570592, 709633,
                                                                       427904, 428688, 583696,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 726328, 0, 3,
                                                                       709633, 571348, 710578,
                                                                       428688, 429472, 584704,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 727588, 0, 3,
                                                                       710578, 572104, 711523,
                                                                       429472, 430256, 585712,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 728848, 0, 3,
                                                                       712468, 573616, 713413,
                                                                       431824, 432608, 586720,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 730108, 0, 3,
                                                                       713413, 574372, 714358,
                                                                       432608, 433392, 587728,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 731368, 0, 3,
                                                                       714358, 575128, 715303,
                                                                       433392, 434176, 588736,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 732628, 0, 3,
                                                                       715303, 575884, 716248,
                                                                       434176, 434960, 589744,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 733888, 0, 3,
                                                                       716248, 576640, 717193,
                                                                       434960, 435744, 590752,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 735148, 0, 3,
                                                                       717193, 577396, 718138,
                                                                       435744, 436528, 591760,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 736408, 0, 3,
                                                                       718138, 578152, 719083,
                                                                       436528, 437312, 592768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 737668, 0, 3,
                                                                       720028, 579664, 721288,
                                                                       438880, 439888, 593776,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 739288, 0, 3,
                                                                       721288, 580672, 722548,
                                                                       439888, 440896, 595072,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 740908, 0, 3,
                                                                       722548, 581680, 723808,
                                                                       440896, 441904, 596368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 742528, 0, 3,
                                                                       723808, 582688, 725068,
                                                                       441904, 442912, 597664,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 744148, 0, 3,
                                                                       725068, 583696, 726328,
                                                                       442912, 443920, 598960,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 745768, 0, 3,
                                                                       726328, 584704, 727588,
                                                                       443920, 444928, 600256,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 747388, 0, 3,
                                                                       728848, 586720, 730108,
                                                                       446944, 447952, 601552,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 749008, 0, 3,
                                                                       730108, 587728, 731368,
                                                                       447952, 448960, 602848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 750628, 0, 3,
                                                                       731368, 588736, 732628,
                                                                       448960, 449968, 604144,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 752248, 0, 3,
                                                                       732628, 589744, 733888,
                                                                       449968, 450976, 605440,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 753868, 0, 3,
                                                                       733888, 590752, 735148,
                                                                       450976, 451984, 606736,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 755488, 0, 3,
                                                                       735148, 591760, 736408,
                                                                       451984, 452992, 608032,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 757108, 0, 3,
                                                                       737668, 593776, 739288,
                                                                       455008, 456268, 609328,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 759133, 0, 3,
                                                                       739288, 595072, 740908,
                                                                       456268, 457528, 610948,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 761158, 0, 3,
                                                                       740908, 596368, 742528,
                                                                       457528, 458788, 612568,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 763183, 0, 3,
                                                                       742528, 597664, 744148,
                                                                       458788, 460048, 614188,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 765208, 0, 3,
                                                                       744148, 598960, 745768,
                                                                       460048, 461308, 615808,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 767233, 0, 3,
                                                                       747388, 601552, 749008,
                                                                       463828, 465088, 617428,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 769258, 0, 3,
                                                                       749008, 602848, 750628,
                                                                       465088, 466348, 619048,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 771283, 0, 3,
                                                                       750628, 604144, 752248,
                                                                       466348, 467608, 620668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 773308, 0, 3,
                                                                       752248, 605440, 753868,
                                                                       467608, 468868, 622288,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 775333, 0, 3,
                                                                       753868, 606736, 755488,
                                                                       468868, 470128, 623908,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 777358, 0, 3,
                                                                       757108, 609328, 759133,
                                                                       472648, 474188, 625528,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 779833, 0, 3,
                                                                       759133, 610948, 761158,
                                                                       474188, 475728, 627508,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 782308, 0, 3,
                                                                       761158, 612568, 763183,
                                                                       475728, 477268, 629488,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 784783, 0, 3,
                                                                       763183, 614188, 765208,
                                                                       477268, 478808, 631468,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 787258, 0, 3,
                                                                       767233, 617428, 769258,
                                                                       481888, 483428, 633448,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 789733, 0, 3,
                                                                       769258, 619048, 771283,
                                                                       483428, 484968, 635428,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 792208, 0, 3,
                                                                       771283, 620668, 773308,
                                                                       484968, 486508, 637408,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 794683, 0, 3,
                                                                       773308, 622288, 775333,
                                                                       486508, 488048, 639388,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 797158, 0, 3,
                                                                       777358, 625528, 779833,
                                                                       491128, 492976, 641368,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 800128, 0, 3,
                                                                       779833, 627508, 782308,
                                                                       492976, 494824, 643744,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 803098, 0, 3,
                                                                       782308, 629488, 784783,
                                                                       494824, 496672, 646120,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 806068, 0, 3,
                                                                       787258, 633448, 789733,
                                                                       500368, 502216, 648496,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 809038, 0, 3,
                                                                       789733, 635428, 792208,
                                                                       502216, 504064, 650872,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 812008, 0, 3,
                                                                       792208, 637408, 794683,
                                                                       504064, 505912, 653248,
                                                                       ncols, gamma, p, q);

                    compute_prim_osl_three_center_electron_repulsion_0(buffer, 814978, 0, 3,
                                                                       797158, 641368, 800128,
                                                                       509608, 511792, 655624,
                                                                       ncols, gamma, p, q);

                    compute_prim_osl_three_center_electron_repulsion_0(buffer, 818488, 0, 3,
                                                                       800128, 643744, 803098,
                                                                       511792, 513976, 658432,
                                                                       ncols, gamma, p, q);

                    compute_prim_osl_three_center_electron_repulsion_0(buffer, 821998, 0, 3,
                                                                       806068, 648496, 809038,
                                                                       518344, 520528, 661240,
                                                                       ncols, gamma, p, q);

                    compute_prim_osl_three_center_electron_repulsion_0(buffer, 825508, 0, 3,
                                                                       809038, 650872, 812008,
                                                                       520528, 522712, 664048,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsl_three_center_electron_repulsion_0(buffer, 829018, 0, 3,
                                                                       814978, 655624, 818488,
                                                                       527080, 529628, 666856,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsl_three_center_electron_repulsion_0(buffer, 833113, 0, 3,
                                                                       821998, 661240, 825508,
                                                                       534724, 537272, 670132,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 837208, 720028, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 838944, 728848, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 840680, 737668, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 842912, 747388, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 845144, 757108, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 847934, 767233, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 850724, 777358, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 854134, 787258, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 857544, 797158, 2970, ncols);

                    simdfunc::contract_primitives(buffer, 861636, 806068, 2970, ncols);

                    simdfunc::contract_primitives(buffer, 865728, 814978, 3510, ncols);

                    simdfunc::contract_primitives(buffer, 870564, 821998, 3510, ncols);

                    simdfunc::contract_primitives(buffer, 875400, 829018, 4095, ncols);

                    simdfunc::contract_primitives(buffer, 881042, 833113, 4095, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 838468, 837208, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 840204, 838944, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 842300, 840680, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 844532, 842912, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 847169, 845144, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 849959, 847934, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 853199, 850724, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 856609, 854134, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 860514, 857544, 66, 1, nmax);

        simdtrf::transform_l_inner(buffer, 864606, 861636, 66, 1, nmax);

        simdtrf::transform_l_inner(buffer, 869238, 865728, 78, 1, nmax);

        simdtrf::transform_l_inner(buffer, 874074, 870564, 78, 1, nmax);

        simdtrf::transform_l_inner(buffer, 879495, 875400, 91, 1, nmax);

        simdtrf::transform_l_inner(buffer, 885137, 881042, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 886684, 838468, 842300, 17,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 888112, 840204, 844532, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 889540, 842300, 847169, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 891376, 844532, 849959, 17,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 893212, 847169, 853199, 17,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 895507, 849959, 856609, 17,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 897802, 853199, 860514, 17,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 900607, 856609, 864606, 17,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 903412, 860514, 869238, 17,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 906778, 864606, 874074, 17,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 910144, 869238, 879495, 17,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 914122, 874074, 885137, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 918100, 886684, 889540, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 920956, 888112, 891376, 17,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 923812, 889540, 893212, 17,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 927484, 891376, 895507, 17,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 931156, 893212, 897802, 17,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 935746, 895507, 900607, 17,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 940336, 897802, 903412, 17,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 945946, 900607, 906778, 17,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 951556, 903412, 910144, 17,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 958288, 906778, 914122, 17,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 965020, 918100, 923812, 17,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 969780, 920956, 927484, 17,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 974540, 923812, 931156, 17,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 980660, 927484, 935746, 17,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 986780, 931156, 940336, 17,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 994430, 935746, 945946, 17,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 1002080, 940336, 951556, 17,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 1011430, 945946, 958288, 17,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 1020780, 965020, 974540, 17,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 1027920, 969780, 980660, 17,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 1035060, 974540, 986780, 17,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 1044240, 980660, 994430, 17,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 1053420, 986780, 1002080, 17,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 1064895, 994430, 1011430, 17,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 1076370, 1020780, 1035060, 17,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 1086366, 1027920, 1044240, 17,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 1096362, 1035060, 1053420, 17,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 1109214, 1044240, 1064895, 17,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 1122066, 1076370, 1096362, 17,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 1135394, 1086366, 1109214, 17,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 1148722, 1135394, 28, 17, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 1148722, 221, nmax);

        simdtrf::transform_i_inner(buffer, 1148722, 1122066, 28, 17, nmax);

        simdtrf::transform_i_outer(values + 2873 * nvalues + n * npairs, nvalues, buffer,
                                   1148722, 221, nmax);
    }

    for (size_t m = 0; m < 5746; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
