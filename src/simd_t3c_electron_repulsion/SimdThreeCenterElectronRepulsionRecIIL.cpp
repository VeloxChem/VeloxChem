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

    const auto nmax = simdfunc::prepare_buffer(buffer, 580552, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 580552, 418607, 23191, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 20,
                                                             ncols, fj, mu, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 7, 8,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 8, 9,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 9, 10,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 10, 11,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 11, 12,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 12, 13,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 124, 0, 3, 13, 14,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 130, 0, 3, 14, 15,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 136, 0, 3, 15, 16,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 16, 17,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 17, 18,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 154, 0, 3, 18, 19,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 160, 0, 3, 19, 20,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 166, 0, 3, 20, 21,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 21, 22,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 22, 23,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 184, 0, 3, 23, 24,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 190, 0, 3, 24, 25,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 196, 0, 3, 25, 26,
                                                                       82, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 28, 31,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 31, 34,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 34, 37,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 37, 40,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 40, 43,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 252, 0, 3, 43, 46,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 262, 0, 3, 46, 49,
                                                                       124, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 49, 52,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 282, 0, 3, 52, 55,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 292, 0, 3, 55, 58,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 302, 0, 3, 58, 61,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 312, 0, 3, 61, 64,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 322, 0, 3, 64, 67,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 332, 0, 3, 67, 70,
                                                                       166, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 342, 0, 3, 70, 73,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 352, 0, 3, 73, 76,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 362, 0, 3, 76, 79,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 372, 0, 3, 79, 82,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 382, 0, 3, 88, 94,
                                                                       202, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 397, 0, 3, 94,
                                                                       100, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 412, 0, 3, 100,
                                                                       106, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 427, 0, 3, 106,
                                                                       112, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 442, 0, 3, 112,
                                                                       118, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 457, 0, 3, 118,
                                                                       124, 252, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 472, 0, 3, 124,
                                                                       130, 262, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 487, 0, 3, 130,
                                                                       136, 272, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 502, 0, 3, 136,
                                                                       142, 282, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 517, 0, 3, 142,
                                                                       148, 292, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 532, 0, 3, 148,
                                                                       154, 302, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 547, 0, 3, 154,
                                                                       160, 312, 322, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 562, 0, 3, 160,
                                                                       166, 322, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 577, 0, 3, 166,
                                                                       172, 332, 342, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 592, 0, 3, 172,
                                                                       178, 342, 352, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 607, 0, 3, 178,
                                                                       184, 352, 362, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 622, 0, 3, 184,
                                                                       190, 362, 372, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 637, 0, 3, 202,
                                                                       212, 382, 397, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 658, 0, 3, 212,
                                                                       222, 397, 412, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 679, 0, 3, 222,
                                                                       232, 412, 427, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 700, 0, 3, 232,
                                                                       242, 427, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 721, 0, 3, 242,
                                                                       252, 442, 457, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 742, 0, 3, 252,
                                                                       262, 457, 472, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 763, 0, 3, 262,
                                                                       272, 472, 487, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 784, 0, 3, 272,
                                                                       282, 487, 502, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 805, 0, 3, 282,
                                                                       292, 502, 517, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 826, 0, 3, 292,
                                                                       302, 517, 532, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 847, 0, 3, 302,
                                                                       312, 532, 547, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 868, 0, 3, 312,
                                                                       322, 547, 562, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 889, 0, 3, 322,
                                                                       332, 562, 577, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 910, 0, 3, 332,
                                                                       342, 577, 592, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 931, 0, 3, 342,
                                                                       352, 592, 607, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 952, 0, 3, 352,
                                                                       362, 607, 622, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 973, 0, 3, 382,
                                                                       397, 637, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1001, 0, 3, 397,
                                                                       412, 658, 679, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1029, 0, 3, 412,
                                                                       427, 679, 700, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1057, 0, 3, 427,
                                                                       442, 700, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1085, 0, 3, 442,
                                                                       457, 721, 742, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1113, 0, 3, 457,
                                                                       472, 742, 763, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1141, 0, 3, 472,
                                                                       487, 763, 784, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1169, 0, 3, 487,
                                                                       502, 784, 805, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1197, 0, 3, 502,
                                                                       517, 805, 826, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1225, 0, 3, 517,
                                                                       532, 826, 847, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1253, 0, 3, 532,
                                                                       547, 847, 868, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1281, 0, 3, 547,
                                                                       562, 868, 889, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1309, 0, 3, 562,
                                                                       577, 889, 910, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1337, 0, 3, 577,
                                                                       592, 910, 931, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1365, 0, 3, 592,
                                                                       607, 931, 952, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1393, 0, 3, 637,
                                                                       658, 973, 1001, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1429, 0, 3, 658,
                                                                       679, 1001, 1029, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1465, 0, 3, 679,
                                                                       700, 1029, 1057, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1501, 0, 3, 700,
                                                                       721, 1057, 1085, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1537, 0, 3, 721,
                                                                       742, 1085, 1113, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1573, 0, 3, 742,
                                                                       763, 1113, 1141, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1609, 0, 3, 763,
                                                                       784, 1141, 1169, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1645, 0, 3, 784,
                                                                       805, 1169, 1197, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1681, 0, 3, 805,
                                                                       826, 1197, 1225, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1717, 0, 3, 826,
                                                                       847, 1225, 1253, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1753, 0, 3, 847,
                                                                       868, 1253, 1281, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1789, 0, 3, 868,
                                                                       889, 1281, 1309, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1825, 0, 3, 889,
                                                                       910, 1309, 1337, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1861, 0, 3, 910,
                                                                       931, 1337, 1365, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1897, 0, 3, 973,
                                                                       1001, 1393, 1429, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1942, 0, 3, 1001,
                                                                       1029, 1429, 1465, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1987, 0, 3, 1029,
                                                                       1057, 1465, 1501, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2032, 0, 3, 1057,
                                                                       1085, 1501, 1537, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2077, 0, 3, 1085,
                                                                       1113, 1537, 1573, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2122, 0, 3, 1113,
                                                                       1141, 1573, 1609, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2167, 0, 3, 1141,
                                                                       1169, 1609, 1645, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2212, 0, 3, 1169,
                                                                       1197, 1645, 1681, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2257, 0, 3, 1197,
                                                                       1225, 1681, 1717, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2302, 0, 3, 1225,
                                                                       1253, 1717, 1753, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2347, 0, 3, 1253,
                                                                       1281, 1753, 1789, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2392, 0, 3, 1281,
                                                                       1309, 1789, 1825, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2437, 0, 3, 1309,
                                                                       1337, 1825, 1861, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2482, 0, 3, 1393,
                                                                       1429, 1897, 1942, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2537, 0, 3, 1429,
                                                                       1465, 1942, 1987, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2592, 0, 3, 1465,
                                                                       1501, 1987, 2032, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2647, 0, 3, 1501,
                                                                       1537, 2032, 2077, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2702, 0, 3, 1537,
                                                                       1573, 2077, 2122, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2757, 0, 3, 1573,
                                                                       1609, 2122, 2167, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2812, 0, 3, 1609,
                                                                       1645, 2167, 2212, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2867, 0, 3, 1645,
                                                                       1681, 2212, 2257, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2922, 0, 3, 1681,
                                                                       1717, 2257, 2302, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2977, 0, 3, 1717,
                                                                       1753, 2302, 2347, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3032, 0, 3, 1753,
                                                                       1789, 2347, 2392, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3087, 0, 3, 1789,
                                                                       1825, 2392, 2437, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3142, 0, 3, 1897,
                                                                       1942, 2482, 2537, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3208, 0, 3, 1942,
                                                                       1987, 2537, 2592, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3274, 0, 3, 1987,
                                                                       2032, 2592, 2647, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3340, 0, 3, 2032,
                                                                       2077, 2647, 2702, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3406, 0, 3, 2077,
                                                                       2122, 2702, 2757, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3472, 0, 3, 2122,
                                                                       2167, 2757, 2812, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3538, 0, 3, 2167,
                                                                       2212, 2812, 2867, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3604, 0, 3, 2212,
                                                                       2257, 2867, 2922, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3670, 0, 3, 2257,
                                                                       2302, 2922, 2977, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3736, 0, 3, 2302,
                                                                       2347, 2977, 3032, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3802, 0, 3, 2347,
                                                                       2392, 3032, 3087, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3868, 0, 3, 2482,
                                                                       2537, 3142, 3208, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3946, 0, 3, 2537,
                                                                       2592, 3208, 3274, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4024, 0, 3, 2592,
                                                                       2647, 3274, 3340, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4102, 0, 3, 2647,
                                                                       2702, 3340, 3406, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4180, 0, 3, 2702,
                                                                       2757, 3406, 3472, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4258, 0, 3, 2757,
                                                                       2812, 3472, 3538, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4336, 0, 3, 2812,
                                                                       2867, 3538, 3604, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4414, 0, 3, 2867,
                                                                       2922, 3604, 3670, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4492, 0, 3, 2922,
                                                                       2977, 3670, 3736, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4570, 0, 3, 2977,
                                                                       3032, 3736, 3802, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4648, 0, 3, 3142,
                                                                       3208, 3868, 3946, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4739, 0, 3, 3208,
                                                                       3274, 3946, 4024, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4830, 0, 3, 3274,
                                                                       3340, 4024, 4102, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4921, 0, 3, 3340,
                                                                       3406, 4102, 4180, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5012, 0, 3, 3406,
                                                                       3472, 4180, 4258, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5103, 0, 3, 3472,
                                                                       3538, 4258, 4336, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5194, 0, 3, 3538,
                                                                       3604, 4336, 4414, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5285, 0, 3, 3604,
                                                                       3670, 4414, 4492, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5376, 0, 3, 3670,
                                                                       3736, 4492, 4570, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5467, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5470, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5473, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5476, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5479, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5482, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5485, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5488, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5491, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5494, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5497, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5500, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5503, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5506, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5509, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5512, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5515, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5518, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5521, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5524, 3, 9, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5533, 3, 10, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5542, 3, 11, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5551, 3, 12, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5560, 3, 13, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5569, 3, 14, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5578, 3, 15, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5587, 3, 16, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5596, 3, 17, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5605, 3, 18, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5614, 3, 19, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5623, 3, 20, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5632, 3, 21, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5641, 3, 22, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5650, 3, 23, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5659, 3, 24, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5668, 3, 25, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5677, 3, 26, 85,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5686, 3, 34, 100,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5704, 3, 37, 106,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5722, 3, 40, 112,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5740, 3, 43, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5758, 3, 46, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5776, 3, 49, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5794, 3, 52, 136,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5812, 3, 55, 142,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5830, 3, 58, 148,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5848, 3, 61, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5866, 3, 64, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5884, 3, 67, 166,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5902, 3, 70, 172,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5920, 3, 73, 178,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5938, 3, 76, 184,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5956, 3, 79, 190,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5974, 3, 82, 196,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5992, 3, 100, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6022, 3, 106, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6052, 3, 112, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6082, 3, 118, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6112, 3, 124, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6142, 3, 130, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6172, 3, 136, 282,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6202, 3, 142, 292,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6232, 3, 148, 302,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6262, 3, 154, 312,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6292, 3, 160, 322,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6322, 3, 166, 332,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6352, 3, 172, 342,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6382, 3, 178, 352,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6412, 3, 184, 362,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6442, 3, 190, 372,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6472, 3, 222, 412,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6517, 3, 232, 427,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6562, 3, 242, 442,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6607, 3, 252, 457,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6652, 3, 262, 472,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6697, 3, 272, 487,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6742, 3, 282, 502,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6787, 3, 292, 517,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6832, 3, 302, 532,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6877, 3, 312, 547,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6922, 3, 322, 562,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6967, 3, 332, 577,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7012, 3, 342, 592,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7057, 3, 352, 607,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7102, 3, 362, 622,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7147, 3, 412, 679,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7210, 3, 427, 700,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7273, 3, 442, 721,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7336, 3, 457, 742,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7399, 3, 472, 763,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7462, 3, 487, 784,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7525, 3, 502, 805,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7588, 3, 517, 826,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7651, 3, 532, 847,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7714, 3, 547, 868,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7777, 3, 562, 889,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7840, 3, 577, 910,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7903, 3, 592, 931,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7966, 3, 607, 952,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8029, 3, 679,
                                                                       1029, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8113, 3, 700,
                                                                       1057, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8197, 3, 721,
                                                                       1085, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8281, 3, 742,
                                                                       1113, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8365, 3, 763,
                                                                       1141, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8449, 3, 784,
                                                                       1169, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8533, 3, 805,
                                                                       1197, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8617, 3, 826,
                                                                       1225, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8701, 3, 847,
                                                                       1253, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8785, 3, 868,
                                                                       1281, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8869, 3, 889,
                                                                       1309, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8953, 3, 910,
                                                                       1337, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9037, 3, 931,
                                                                       1365, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9121, 3, 1029,
                                                                       1465, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9229, 3, 1057,
                                                                       1501, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9337, 3, 1085,
                                                                       1537, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9445, 3, 1113,
                                                                       1573, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9553, 3, 1141,
                                                                       1609, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9661, 3, 1169,
                                                                       1645, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9769, 3, 1197,
                                                                       1681, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9877, 3, 1225,
                                                                       1717, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9985, 3, 1253,
                                                                       1753, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10093, 3, 1281,
                                                                       1789, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10201, 3, 1309,
                                                                       1825, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10309, 3, 1337,
                                                                       1861, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10417, 3, 1465,
                                                                       1987, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10552, 3, 1501,
                                                                       2032, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10687, 3, 1537,
                                                                       2077, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10822, 3, 1573,
                                                                       2122, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10957, 3, 1609,
                                                                       2167, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11092, 3, 1645,
                                                                       2212, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11227, 3, 1681,
                                                                       2257, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11362, 3, 1717,
                                                                       2302, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11497, 3, 1753,
                                                                       2347, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11632, 3, 1789,
                                                                       2392, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11767, 3, 1825,
                                                                       2437, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11902, 3, 1987,
                                                                       2592, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12067, 3, 2032,
                                                                       2647, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12232, 3, 2077,
                                                                       2702, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12397, 3, 2122,
                                                                       2757, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12562, 3, 2167,
                                                                       2812, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12727, 3, 2212,
                                                                       2867, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12892, 3, 2257,
                                                                       2922, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13057, 3, 2302,
                                                                       2977, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13222, 3, 2347,
                                                                       3032, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13387, 3, 2392,
                                                                       3087, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13552, 3, 2592,
                                                                       3274, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13750, 3, 2647,
                                                                       3340, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13948, 3, 2702,
                                                                       3406, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14146, 3, 2757,
                                                                       3472, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14344, 3, 2812,
                                                                       3538, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14542, 3, 2867,
                                                                       3604, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14740, 3, 2922,
                                                                       3670, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14938, 3, 2977,
                                                                       3736, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15136, 3, 3032,
                                                                       3802, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15334, 3, 3274,
                                                                       4024, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15568, 3, 3340,
                                                                       4102, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15802, 3, 3406,
                                                                       4180, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16036, 3, 3472,
                                                                       4258, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16270, 3, 3538,
                                                                       4336, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16504, 3, 3604,
                                                                       4414, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16738, 3, 3670,
                                                                       4492, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16972, 3, 3736,
                                                                       4570, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17206, 3, 4024,
                                                                       4830, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17479, 3, 4102,
                                                                       4921, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17752, 3, 4180,
                                                                       5012, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 18025, 3, 4258,
                                                                       5103, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 18298, 3, 4336,
                                                                       5194, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 18571, 3, 4414,
                                                                       5285, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 18844, 3, 4492,
                                                                       5376, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19117, 3, 7, 8,
                                                                       5467, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19123, 3, 8, 9,
                                                                       5470, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19129, 3, 9, 10,
                                                                       5473, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19135, 3, 10, 11,
                                                                       5476, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19141, 3, 11, 12,
                                                                       5479, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19147, 3, 12, 13,
                                                                       5482, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19153, 3, 13, 14,
                                                                       5485, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19159, 3, 14, 15,
                                                                       5488, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19165, 3, 15, 16,
                                                                       5491, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19171, 3, 16, 17,
                                                                       5494, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19177, 3, 17, 18,
                                                                       5497, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19183, 3, 18, 19,
                                                                       5500, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19189, 3, 19, 20,
                                                                       5503, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19195, 3, 20, 21,
                                                                       5506, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19201, 3, 21, 22,
                                                                       5509, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19207, 3, 22, 23,
                                                                       5512, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19213, 3, 23, 24,
                                                                       5515, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19219, 3, 24, 25,
                                                                       5518, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19225, 3, 25, 26,
                                                                       5521, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19231, 0, 3,
                                                                       19117, 5467, 19123, 5524,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19249, 0, 3,
                                                                       19123, 5470, 19129, 5533,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19267, 0, 3,
                                                                       19129, 5473, 19135, 5542,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19285, 0, 3,
                                                                       19135, 5476, 19141, 5551,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19303, 0, 3,
                                                                       19141, 5479, 19147, 5560,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19321, 0, 3,
                                                                       19147, 5482, 19153, 5569,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19339, 0, 3,
                                                                       19153, 5485, 19159, 5578,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19357, 0, 3,
                                                                       19159, 5488, 19165, 5587,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19375, 0, 3,
                                                                       19165, 5491, 19171, 5596,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19393, 0, 3,
                                                                       19171, 5494, 19177, 5605,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19411, 0, 3,
                                                                       19177, 5497, 19183, 5614,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19429, 0, 3,
                                                                       19183, 5500, 19189, 5623,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19447, 0, 3,
                                                                       19189, 5503, 19195, 5632,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19465, 0, 3,
                                                                       19195, 5506, 19201, 5641,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19483, 0, 3,
                                                                       19201, 5509, 19207, 5650,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19501, 0, 3,
                                                                       19207, 5512, 19213, 5659,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19519, 0, 3,
                                                                       19213, 5515, 19219, 5668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19537, 0, 3,
                                                                       19219, 5518, 19225, 5677,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19555, 0, 3,
                                                                       19231, 5524, 19249, 88,
                                                                       94, 5686, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19591, 0, 3,
                                                                       19249, 5533, 19267, 94,
                                                                       100, 5704, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19627, 0, 3,
                                                                       19267, 5542, 19285, 100,
                                                                       106, 5722, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19663, 0, 3,
                                                                       19285, 5551, 19303, 106,
                                                                       112, 5740, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19699, 0, 3,
                                                                       19303, 5560, 19321, 112,
                                                                       118, 5758, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19735, 0, 3,
                                                                       19321, 5569, 19339, 118,
                                                                       124, 5776, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19771, 0, 3,
                                                                       19339, 5578, 19357, 124,
                                                                       130, 5794, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19807, 0, 3,
                                                                       19357, 5587, 19375, 130,
                                                                       136, 5812, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19843, 0, 3,
                                                                       19375, 5596, 19393, 136,
                                                                       142, 5830, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19879, 0, 3,
                                                                       19393, 5605, 19411, 142,
                                                                       148, 5848, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19915, 0, 3,
                                                                       19411, 5614, 19429, 148,
                                                                       154, 5866, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19951, 0, 3,
                                                                       19429, 5623, 19447, 154,
                                                                       160, 5884, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19987, 0, 3,
                                                                       19447, 5632, 19465, 160,
                                                                       166, 5902, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 20023, 0, 3,
                                                                       19465, 5641, 19483, 166,
                                                                       172, 5920, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 20059, 0, 3,
                                                                       19483, 5650, 19501, 172,
                                                                       178, 5938, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 20095, 0, 3,
                                                                       19501, 5659, 19519, 178,
                                                                       184, 5956, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 20131, 0, 3,
                                                                       19519, 5668, 19537, 184,
                                                                       190, 5974, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20167, 0, 3,
                                                                       19555, 5686, 19591, 202,
                                                                       212, 5992, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20227, 0, 3,
                                                                       19591, 5704, 19627, 212,
                                                                       222, 6022, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20287, 0, 3,
                                                                       19627, 5722, 19663, 222,
                                                                       232, 6052, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20347, 0, 3,
                                                                       19663, 5740, 19699, 232,
                                                                       242, 6082, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20407, 0, 3,
                                                                       19699, 5758, 19735, 242,
                                                                       252, 6112, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20467, 0, 3,
                                                                       19735, 5776, 19771, 252,
                                                                       262, 6142, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20527, 0, 3,
                                                                       19771, 5794, 19807, 262,
                                                                       272, 6172, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20587, 0, 3,
                                                                       19807, 5812, 19843, 272,
                                                                       282, 6202, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20647, 0, 3,
                                                                       19843, 5830, 19879, 282,
                                                                       292, 6232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20707, 0, 3,
                                                                       19879, 5848, 19915, 292,
                                                                       302, 6262, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20767, 0, 3,
                                                                       19915, 5866, 19951, 302,
                                                                       312, 6292, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20827, 0, 3,
                                                                       19951, 5884, 19987, 312,
                                                                       322, 6322, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20887, 0, 3,
                                                                       19987, 5902, 20023, 322,
                                                                       332, 6352, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20947, 0, 3,
                                                                       20023, 5920, 20059, 332,
                                                                       342, 6382, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 21007, 0, 3,
                                                                       20059, 5938, 20095, 342,
                                                                       352, 6412, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 21067, 0, 3,
                                                                       20095, 5956, 20131, 352,
                                                                       362, 6442, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21127, 0, 3,
                                                                       20167, 5992, 20227, 382,
                                                                       397, 6472, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21217, 0, 3,
                                                                       20227, 6022, 20287, 397,
                                                                       412, 6517, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21307, 0, 3,
                                                                       20287, 6052, 20347, 412,
                                                                       427, 6562, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21397, 0, 3,
                                                                       20347, 6082, 20407, 427,
                                                                       442, 6607, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21487, 0, 3,
                                                                       20407, 6112, 20467, 442,
                                                                       457, 6652, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21577, 0, 3,
                                                                       20467, 6142, 20527, 457,
                                                                       472, 6697, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21667, 0, 3,
                                                                       20527, 6172, 20587, 472,
                                                                       487, 6742, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21757, 0, 3,
                                                                       20587, 6202, 20647, 487,
                                                                       502, 6787, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21847, 0, 3,
                                                                       20647, 6232, 20707, 502,
                                                                       517, 6832, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21937, 0, 3,
                                                                       20707, 6262, 20767, 517,
                                                                       532, 6877, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22027, 0, 3,
                                                                       20767, 6292, 20827, 532,
                                                                       547, 6922, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22117, 0, 3,
                                                                       20827, 6322, 20887, 547,
                                                                       562, 6967, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22207, 0, 3,
                                                                       20887, 6352, 20947, 562,
                                                                       577, 7012, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22297, 0, 3,
                                                                       20947, 6382, 21007, 577,
                                                                       592, 7057, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22387, 0, 3,
                                                                       21007, 6412, 21067, 592,
                                                                       607, 7102, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22477, 0, 3,
                                                                       21127, 6472, 21217, 637,
                                                                       658, 7147, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22603, 0, 3,
                                                                       21217, 6517, 21307, 658,
                                                                       679, 7210, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22729, 0, 3,
                                                                       21307, 6562, 21397, 679,
                                                                       700, 7273, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22855, 0, 3,
                                                                       21397, 6607, 21487, 700,
                                                                       721, 7336, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22981, 0, 3,
                                                                       21487, 6652, 21577, 721,
                                                                       742, 7399, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23107, 0, 3,
                                                                       21577, 6697, 21667, 742,
                                                                       763, 7462, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23233, 0, 3,
                                                                       21667, 6742, 21757, 763,
                                                                       784, 7525, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23359, 0, 3,
                                                                       21757, 6787, 21847, 784,
                                                                       805, 7588, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23485, 0, 3,
                                                                       21847, 6832, 21937, 805,
                                                                       826, 7651, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23611, 0, 3,
                                                                       21937, 6877, 22027, 826,
                                                                       847, 7714, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23737, 0, 3,
                                                                       22027, 6922, 22117, 847,
                                                                       868, 7777, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23863, 0, 3,
                                                                       22117, 6967, 22207, 868,
                                                                       889, 7840, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23989, 0, 3,
                                                                       22207, 7012, 22297, 889,
                                                                       910, 7903, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24115, 0, 3,
                                                                       22297, 7057, 22387, 910,
                                                                       931, 7966, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24241, 0, 3,
                                                                       22477, 7147, 22603, 973,
                                                                       1001, 8029, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24409, 0, 3,
                                                                       22603, 7210, 22729, 1001,
                                                                       1029, 8113, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24577, 0, 3,
                                                                       22729, 7273, 22855, 1029,
                                                                       1057, 8197, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24745, 0, 3,
                                                                       22855, 7336, 22981, 1057,
                                                                       1085, 8281, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24913, 0, 3,
                                                                       22981, 7399, 23107, 1085,
                                                                       1113, 8365, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25081, 0, 3,
                                                                       23107, 7462, 23233, 1113,
                                                                       1141, 8449, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25249, 0, 3,
                                                                       23233, 7525, 23359, 1141,
                                                                       1169, 8533, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25417, 0, 3,
                                                                       23359, 7588, 23485, 1169,
                                                                       1197, 8617, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25585, 0, 3,
                                                                       23485, 7651, 23611, 1197,
                                                                       1225, 8701, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25753, 0, 3,
                                                                       23611, 7714, 23737, 1225,
                                                                       1253, 8785, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25921, 0, 3,
                                                                       23737, 7777, 23863, 1253,
                                                                       1281, 8869, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26089, 0, 3,
                                                                       23863, 7840, 23989, 1281,
                                                                       1309, 8953, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26257, 0, 3,
                                                                       23989, 7903, 24115, 1309,
                                                                       1337, 9037, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26425, 0, 3,
                                                                       24241, 8029, 24409, 1393,
                                                                       1429, 9121, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26641, 0, 3,
                                                                       24409, 8113, 24577, 1429,
                                                                       1465, 9229, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26857, 0, 3,
                                                                       24577, 8197, 24745, 1465,
                                                                       1501, 9337, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27073, 0, 3,
                                                                       24745, 8281, 24913, 1501,
                                                                       1537, 9445, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27289, 0, 3,
                                                                       24913, 8365, 25081, 1537,
                                                                       1573, 9553, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27505, 0, 3,
                                                                       25081, 8449, 25249, 1573,
                                                                       1609, 9661, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27721, 0, 3,
                                                                       25249, 8533, 25417, 1609,
                                                                       1645, 9769, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27937, 0, 3,
                                                                       25417, 8617, 25585, 1645,
                                                                       1681, 9877, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28153, 0, 3,
                                                                       25585, 8701, 25753, 1681,
                                                                       1717, 9985, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28369, 0, 3,
                                                                       25753, 8785, 25921, 1717,
                                                                       1753, 10093, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28585, 0, 3,
                                                                       25921, 8869, 26089, 1753,
                                                                       1789, 10201, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28801, 0, 3,
                                                                       26089, 8953, 26257, 1789,
                                                                       1825, 10309, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29017, 0, 3,
                                                                       26425, 9121, 26641, 1897,
                                                                       1942, 10417, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29287, 0, 3,
                                                                       26641, 9229, 26857, 1942,
                                                                       1987, 10552, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29557, 0, 3,
                                                                       26857, 9337, 27073, 1987,
                                                                       2032, 10687, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29827, 0, 3,
                                                                       27073, 9445, 27289, 2032,
                                                                       2077, 10822, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30097, 0, 3,
                                                                       27289, 9553, 27505, 2077,
                                                                       2122, 10957, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30367, 0, 3,
                                                                       27505, 9661, 27721, 2122,
                                                                       2167, 11092, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30637, 0, 3,
                                                                       27721, 9769, 27937, 2167,
                                                                       2212, 11227, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30907, 0, 3,
                                                                       27937, 9877, 28153, 2212,
                                                                       2257, 11362, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31177, 0, 3,
                                                                       28153, 9985, 28369, 2257,
                                                                       2302, 11497, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31447, 0, 3,
                                                                       28369, 10093, 28585, 2302,
                                                                       2347, 11632, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31717, 0, 3,
                                                                       28585, 10201, 28801, 2347,
                                                                       2392, 11767, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31987, 0, 3,
                                                                       29017, 10417, 29287, 2482,
                                                                       2537, 11902, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32317, 0, 3,
                                                                       29287, 10552, 29557, 2537,
                                                                       2592, 12067, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32647, 0, 3,
                                                                       29557, 10687, 29827, 2592,
                                                                       2647, 12232, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32977, 0, 3,
                                                                       29827, 10822, 30097, 2647,
                                                                       2702, 12397, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 33307, 0, 3,
                                                                       30097, 10957, 30367, 2702,
                                                                       2757, 12562, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 33637, 0, 3,
                                                                       30367, 11092, 30637, 2757,
                                                                       2812, 12727, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 33967, 0, 3,
                                                                       30637, 11227, 30907, 2812,
                                                                       2867, 12892, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34297, 0, 3,
                                                                       30907, 11362, 31177, 2867,
                                                                       2922, 13057, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34627, 0, 3,
                                                                       31177, 11497, 31447, 2922,
                                                                       2977, 13222, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34957, 0, 3,
                                                                       31447, 11632, 31717, 2977,
                                                                       3032, 13387, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 35287, 0, 3,
                                                                       31987, 11902, 32317, 3142,
                                                                       3208, 13552, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 35683, 0, 3,
                                                                       32317, 12067, 32647, 3208,
                                                                       3274, 13750, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 36079, 0, 3,
                                                                       32647, 12232, 32977, 3274,
                                                                       3340, 13948, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 36475, 0, 3,
                                                                       32977, 12397, 33307, 3340,
                                                                       3406, 14146, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 36871, 0, 3,
                                                                       33307, 12562, 33637, 3406,
                                                                       3472, 14344, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 37267, 0, 3,
                                                                       33637, 12727, 33967, 3472,
                                                                       3538, 14542, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 37663, 0, 3,
                                                                       33967, 12892, 34297, 3538,
                                                                       3604, 14740, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 38059, 0, 3,
                                                                       34297, 13057, 34627, 3604,
                                                                       3670, 14938, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 38455, 0, 3,
                                                                       34627, 13222, 34957, 3670,
                                                                       3736, 15136, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 38851, 0, 3,
                                                                       35287, 13552, 35683, 3868,
                                                                       3946, 15334, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 39319, 0, 3,
                                                                       35683, 13750, 36079, 3946,
                                                                       4024, 15568, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 39787, 0, 3,
                                                                       36079, 13948, 36475, 4024,
                                                                       4102, 15802, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 40255, 0, 3,
                                                                       36475, 14146, 36871, 4102,
                                                                       4180, 16036, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 40723, 0, 3,
                                                                       36871, 14344, 37267, 4180,
                                                                       4258, 16270, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 41191, 0, 3,
                                                                       37267, 14542, 37663, 4258,
                                                                       4336, 16504, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 41659, 0, 3,
                                                                       37663, 14740, 38059, 4336,
                                                                       4414, 16738, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 42127, 0, 3,
                                                                       38059, 14938, 38455, 4414,
                                                                       4492, 16972, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 42595, 0, 3,
                                                                       38851, 15334, 39319, 4648,
                                                                       4739, 17206, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 43141, 0, 3,
                                                                       39319, 15568, 39787, 4739,
                                                                       4830, 17479, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 43687, 0, 3,
                                                                       39787, 15802, 40255, 4830,
                                                                       4921, 17752, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 44233, 0, 3,
                                                                       40255, 16036, 40723, 4921,
                                                                       5012, 18025, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 44779, 0, 3,
                                                                       40723, 16270, 41191, 5012,
                                                                       5103, 18298, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 45325, 0, 3,
                                                                       41191, 16504, 41659, 5103,
                                                                       5194, 18571, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 45871, 0, 3,
                                                                       41659, 16738, 42127, 5194,
                                                                       5285, 18844, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46417, 3, 5467,
                                                                       5470, 19129, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46427, 3, 5470,
                                                                       5473, 19135, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46437, 3, 5473,
                                                                       5476, 19141, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46447, 3, 5476,
                                                                       5479, 19147, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46457, 3, 5479,
                                                                       5482, 19153, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46467, 3, 5482,
                                                                       5485, 19159, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46477, 3, 5485,
                                                                       5488, 19165, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46487, 3, 5488,
                                                                       5491, 19171, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46497, 3, 5491,
                                                                       5494, 19177, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46507, 3, 5494,
                                                                       5497, 19183, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46517, 3, 5497,
                                                                       5500, 19189, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46527, 3, 5500,
                                                                       5503, 19195, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46537, 3, 5503,
                                                                       5506, 19201, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46547, 3, 5506,
                                                                       5509, 19207, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46557, 3, 5509,
                                                                       5512, 19213, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46567, 3, 5512,
                                                                       5515, 19219, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46577, 3, 5515,
                                                                       5518, 19225, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46587, 0, 3,
                                                                       46417, 19129, 46427, 5524,
                                                                       5533, 19267, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46617, 0, 3,
                                                                       46427, 19135, 46437, 5533,
                                                                       5542, 19285, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46647, 0, 3,
                                                                       46437, 19141, 46447, 5542,
                                                                       5551, 19303, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46677, 0, 3,
                                                                       46447, 19147, 46457, 5551,
                                                                       5560, 19321, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46707, 0, 3,
                                                                       46457, 19153, 46467, 5560,
                                                                       5569, 19339, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46737, 0, 3,
                                                                       46467, 19159, 46477, 5569,
                                                                       5578, 19357, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46767, 0, 3,
                                                                       46477, 19165, 46487, 5578,
                                                                       5587, 19375, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46797, 0, 3,
                                                                       46487, 19171, 46497, 5587,
                                                                       5596, 19393, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46827, 0, 3,
                                                                       46497, 19177, 46507, 5596,
                                                                       5605, 19411, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46857, 0, 3,
                                                                       46507, 19183, 46517, 5605,
                                                                       5614, 19429, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46887, 0, 3,
                                                                       46517, 19189, 46527, 5614,
                                                                       5623, 19447, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46917, 0, 3,
                                                                       46527, 19195, 46537, 5623,
                                                                       5632, 19465, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46947, 0, 3,
                                                                       46537, 19201, 46547, 5632,
                                                                       5641, 19483, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46977, 0, 3,
                                                                       46547, 19207, 46557, 5641,
                                                                       5650, 19501, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47007, 0, 3,
                                                                       46557, 19213, 46567, 5650,
                                                                       5659, 19519, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47037, 0, 3,
                                                                       46567, 19219, 46577, 5659,
                                                                       5668, 19537, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47067, 0, 3,
                                                                       46587, 19267, 46617, 5686,
                                                                       5704, 19627, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47127, 0, 3,
                                                                       46617, 19285, 46647, 5704,
                                                                       5722, 19663, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47187, 0, 3,
                                                                       46647, 19303, 46677, 5722,
                                                                       5740, 19699, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47247, 0, 3,
                                                                       46677, 19321, 46707, 5740,
                                                                       5758, 19735, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47307, 0, 3,
                                                                       46707, 19339, 46737, 5758,
                                                                       5776, 19771, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47367, 0, 3,
                                                                       46737, 19357, 46767, 5776,
                                                                       5794, 19807, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47427, 0, 3,
                                                                       46767, 19375, 46797, 5794,
                                                                       5812, 19843, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47487, 0, 3,
                                                                       46797, 19393, 46827, 5812,
                                                                       5830, 19879, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47547, 0, 3,
                                                                       46827, 19411, 46857, 5830,
                                                                       5848, 19915, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47607, 0, 3,
                                                                       46857, 19429, 46887, 5848,
                                                                       5866, 19951, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47667, 0, 3,
                                                                       46887, 19447, 46917, 5866,
                                                                       5884, 19987, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47727, 0, 3,
                                                                       46917, 19465, 46947, 5884,
                                                                       5902, 20023, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47787, 0, 3,
                                                                       46947, 19483, 46977, 5902,
                                                                       5920, 20059, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47847, 0, 3,
                                                                       46977, 19501, 47007, 5920,
                                                                       5938, 20095, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47907, 0, 3,
                                                                       47007, 19519, 47037, 5938,
                                                                       5956, 20131, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47967, 0, 3,
                                                                       47067, 19627, 47127, 5992,
                                                                       6022, 20287, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48067, 0, 3,
                                                                       47127, 19663, 47187, 6022,
                                                                       6052, 20347, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48167, 0, 3,
                                                                       47187, 19699, 47247, 6052,
                                                                       6082, 20407, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48267, 0, 3,
                                                                       47247, 19735, 47307, 6082,
                                                                       6112, 20467, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48367, 0, 3,
                                                                       47307, 19771, 47367, 6112,
                                                                       6142, 20527, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48467, 0, 3,
                                                                       47367, 19807, 47427, 6142,
                                                                       6172, 20587, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48567, 0, 3,
                                                                       47427, 19843, 47487, 6172,
                                                                       6202, 20647, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48667, 0, 3,
                                                                       47487, 19879, 47547, 6202,
                                                                       6232, 20707, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48767, 0, 3,
                                                                       47547, 19915, 47607, 6232,
                                                                       6262, 20767, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48867, 0, 3,
                                                                       47607, 19951, 47667, 6262,
                                                                       6292, 20827, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48967, 0, 3,
                                                                       47667, 19987, 47727, 6292,
                                                                       6322, 20887, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 49067, 0, 3,
                                                                       47727, 20023, 47787, 6322,
                                                                       6352, 20947, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 49167, 0, 3,
                                                                       47787, 20059, 47847, 6352,
                                                                       6382, 21007, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 49267, 0, 3,
                                                                       47847, 20095, 47907, 6382,
                                                                       6412, 21067, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49367, 0, 3,
                                                                       47967, 20287, 48067, 6472,
                                                                       6517, 21307, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49517, 0, 3,
                                                                       48067, 20347, 48167, 6517,
                                                                       6562, 21397, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49667, 0, 3,
                                                                       48167, 20407, 48267, 6562,
                                                                       6607, 21487, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49817, 0, 3,
                                                                       48267, 20467, 48367, 6607,
                                                                       6652, 21577, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49967, 0, 3,
                                                                       48367, 20527, 48467, 6652,
                                                                       6697, 21667, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50117, 0, 3,
                                                                       48467, 20587, 48567, 6697,
                                                                       6742, 21757, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50267, 0, 3,
                                                                       48567, 20647, 48667, 6742,
                                                                       6787, 21847, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50417, 0, 3,
                                                                       48667, 20707, 48767, 6787,
                                                                       6832, 21937, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50567, 0, 3,
                                                                       48767, 20767, 48867, 6832,
                                                                       6877, 22027, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50717, 0, 3,
                                                                       48867, 20827, 48967, 6877,
                                                                       6922, 22117, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50867, 0, 3,
                                                                       48967, 20887, 49067, 6922,
                                                                       6967, 22207, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 51017, 0, 3,
                                                                       49067, 20947, 49167, 6967,
                                                                       7012, 22297, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 51167, 0, 3,
                                                                       49167, 21007, 49267, 7012,
                                                                       7057, 22387, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51317, 0, 3,
                                                                       49367, 21307, 49517, 7147,
                                                                       7210, 22729, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51527, 0, 3,
                                                                       49517, 21397, 49667, 7210,
                                                                       7273, 22855, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51737, 0, 3,
                                                                       49667, 21487, 49817, 7273,
                                                                       7336, 22981, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51947, 0, 3,
                                                                       49817, 21577, 49967, 7336,
                                                                       7399, 23107, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52157, 0, 3,
                                                                       49967, 21667, 50117, 7399,
                                                                       7462, 23233, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52367, 0, 3,
                                                                       50117, 21757, 50267, 7462,
                                                                       7525, 23359, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52577, 0, 3,
                                                                       50267, 21847, 50417, 7525,
                                                                       7588, 23485, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52787, 0, 3,
                                                                       50417, 21937, 50567, 7588,
                                                                       7651, 23611, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52997, 0, 3,
                                                                       50567, 22027, 50717, 7651,
                                                                       7714, 23737, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53207, 0, 3,
                                                                       50717, 22117, 50867, 7714,
                                                                       7777, 23863, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53417, 0, 3,
                                                                       50867, 22207, 51017, 7777,
                                                                       7840, 23989, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53627, 0, 3,
                                                                       51017, 22297, 51167, 7840,
                                                                       7903, 24115, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 53837, 0, 3,
                                                                       51317, 22729, 51527, 8029,
                                                                       8113, 24577, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54117, 0, 3,
                                                                       51527, 22855, 51737, 8113,
                                                                       8197, 24745, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54397, 0, 3,
                                                                       51737, 22981, 51947, 8197,
                                                                       8281, 24913, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54677, 0, 3,
                                                                       51947, 23107, 52157, 8281,
                                                                       8365, 25081, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54957, 0, 3,
                                                                       52157, 23233, 52367, 8365,
                                                                       8449, 25249, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55237, 0, 3,
                                                                       52367, 23359, 52577, 8449,
                                                                       8533, 25417, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55517, 0, 3,
                                                                       52577, 23485, 52787, 8533,
                                                                       8617, 25585, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55797, 0, 3,
                                                                       52787, 23611, 52997, 8617,
                                                                       8701, 25753, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56077, 0, 3,
                                                                       52997, 23737, 53207, 8701,
                                                                       8785, 25921, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56357, 0, 3,
                                                                       53207, 23863, 53417, 8785,
                                                                       8869, 26089, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56637, 0, 3,
                                                                       53417, 23989, 53627, 8869,
                                                                       8953, 26257, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 56917, 0, 3,
                                                                       53837, 24577, 54117, 9121,
                                                                       9229, 26857, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57277, 0, 3,
                                                                       54117, 24745, 54397, 9229,
                                                                       9337, 27073, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57637, 0, 3,
                                                                       54397, 24913, 54677, 9337,
                                                                       9445, 27289, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57997, 0, 3,
                                                                       54677, 25081, 54957, 9445,
                                                                       9553, 27505, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58357, 0, 3,
                                                                       54957, 25249, 55237, 9553,
                                                                       9661, 27721, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58717, 0, 3,
                                                                       55237, 25417, 55517, 9661,
                                                                       9769, 27937, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59077, 0, 3,
                                                                       55517, 25585, 55797, 9769,
                                                                       9877, 28153, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59437, 0, 3,
                                                                       55797, 25753, 56077, 9877,
                                                                       9985, 28369, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59797, 0, 3,
                                                                       56077, 25921, 56357, 9985,
                                                                       10093, 28585, ncols,
                                                                       gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 60157, 0, 3,
                                                                       56357, 26089, 56637,
                                                                       10093, 10201, 28801,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 60517, 0, 3,
                                                                       56917, 26857, 57277,
                                                                       10417, 10552, 29557,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 60967, 0, 3,
                                                                       57277, 27073, 57637,
                                                                       10552, 10687, 29827,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 61417, 0, 3,
                                                                       57637, 27289, 57997,
                                                                       10687, 10822, 30097,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 61867, 0, 3,
                                                                       57997, 27505, 58357,
                                                                       10822, 10957, 30367,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 62317, 0, 3,
                                                                       58357, 27721, 58717,
                                                                       10957, 11092, 30637,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 62767, 0, 3,
                                                                       58717, 27937, 59077,
                                                                       11092, 11227, 30907,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 63217, 0, 3,
                                                                       59077, 28153, 59437,
                                                                       11227, 11362, 31177,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 63667, 0, 3,
                                                                       59437, 28369, 59797,
                                                                       11362, 11497, 31447,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 64117, 0, 3,
                                                                       59797, 28585, 60157,
                                                                       11497, 11632, 31717,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 64567, 0, 3,
                                                                       60517, 29557, 60967,
                                                                       11902, 12067, 32647,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 65117, 0, 3,
                                                                       60967, 29827, 61417,
                                                                       12067, 12232, 32977,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 65667, 0, 3,
                                                                       61417, 30097, 61867,
                                                                       12232, 12397, 33307,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 66217, 0, 3,
                                                                       61867, 30367, 62317,
                                                                       12397, 12562, 33637,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 66767, 0, 3,
                                                                       62317, 30637, 62767,
                                                                       12562, 12727, 33967,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 67317, 0, 3,
                                                                       62767, 30907, 63217,
                                                                       12727, 12892, 34297,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 67867, 0, 3,
                                                                       63217, 31177, 63667,
                                                                       12892, 13057, 34627,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 68417, 0, 3,
                                                                       63667, 31447, 64117,
                                                                       13057, 13222, 34957,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 68967, 0, 3,
                                                                       64567, 32647, 65117,
                                                                       13552, 13750, 36079,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 69627, 0, 3,
                                                                       65117, 32977, 65667,
                                                                       13750, 13948, 36475,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 70287, 0, 3,
                                                                       65667, 33307, 66217,
                                                                       13948, 14146, 36871,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 70947, 0, 3,
                                                                       66217, 33637, 66767,
                                                                       14146, 14344, 37267,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 71607, 0, 3,
                                                                       66767, 33967, 67317,
                                                                       14344, 14542, 37663,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 72267, 0, 3,
                                                                       67317, 34297, 67867,
                                                                       14542, 14740, 38059,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 72927, 0, 3,
                                                                       67867, 34627, 68417,
                                                                       14740, 14938, 38455,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 73587, 0, 3,
                                                                       68967, 36079, 69627,
                                                                       15334, 15568, 39787,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 74367, 0, 3,
                                                                       69627, 36475, 70287,
                                                                       15568, 15802, 40255,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 75147, 0, 3,
                                                                       70287, 36871, 70947,
                                                                       15802, 16036, 40723,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 75927, 0, 3,
                                                                       70947, 37267, 71607,
                                                                       16036, 16270, 41191,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 76707, 0, 3,
                                                                       71607, 37663, 72267,
                                                                       16270, 16504, 41659,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 77487, 0, 3,
                                                                       72267, 38059, 72927,
                                                                       16504, 16738, 42127,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 78267, 0, 3,
                                                                       73587, 39787, 74367,
                                                                       17206, 17479, 43687,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 79177, 0, 3,
                                                                       74367, 40255, 75147,
                                                                       17479, 17752, 44233,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 80087, 0, 3,
                                                                       75147, 40723, 75927,
                                                                       17752, 18025, 44779,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 80997, 0, 3,
                                                                       75927, 41191, 76707,
                                                                       18025, 18298, 45325,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 81907, 0, 3,
                                                                       76707, 41659, 77487,
                                                                       18298, 18571, 45871,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82817, 3, 19117,
                                                                       19123, 46417, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82832, 3, 19123,
                                                                       19129, 46427, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82847, 3, 19129,
                                                                       19135, 46437, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82862, 3, 19135,
                                                                       19141, 46447, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82877, 3, 19141,
                                                                       19147, 46457, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82892, 3, 19147,
                                                                       19153, 46467, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82907, 3, 19153,
                                                                       19159, 46477, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82922, 3, 19159,
                                                                       19165, 46487, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82937, 3, 19165,
                                                                       19171, 46497, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82952, 3, 19171,
                                                                       19177, 46507, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82967, 3, 19177,
                                                                       19183, 46517, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82982, 3, 19183,
                                                                       19189, 46527, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82997, 3, 19189,
                                                                       19195, 46537, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 83012, 3, 19195,
                                                                       19201, 46547, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 83027, 3, 19201,
                                                                       19207, 46557, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 83042, 3, 19207,
                                                                       19213, 46567, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 83057, 3, 19213,
                                                                       19219, 46577, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83072, 0, 3,
                                                                       82817, 46417, 82832,
                                                                       19231, 19249, 46587,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83117, 0, 3,
                                                                       82832, 46427, 82847,
                                                                       19249, 19267, 46617,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83162, 0, 3,
                                                                       82847, 46437, 82862,
                                                                       19267, 19285, 46647,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83207, 0, 3,
                                                                       82862, 46447, 82877,
                                                                       19285, 19303, 46677,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83252, 0, 3,
                                                                       82877, 46457, 82892,
                                                                       19303, 19321, 46707,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83297, 0, 3,
                                                                       82892, 46467, 82907,
                                                                       19321, 19339, 46737,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83342, 0, 3,
                                                                       82907, 46477, 82922,
                                                                       19339, 19357, 46767,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83387, 0, 3,
                                                                       82922, 46487, 82937,
                                                                       19357, 19375, 46797,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83432, 0, 3,
                                                                       82937, 46497, 82952,
                                                                       19375, 19393, 46827,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83477, 0, 3,
                                                                       82952, 46507, 82967,
                                                                       19393, 19411, 46857,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83522, 0, 3,
                                                                       82967, 46517, 82982,
                                                                       19411, 19429, 46887,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83567, 0, 3,
                                                                       82982, 46527, 82997,
                                                                       19429, 19447, 46917,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83612, 0, 3,
                                                                       82997, 46537, 83012,
                                                                       19447, 19465, 46947,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83657, 0, 3,
                                                                       83012, 46547, 83027,
                                                                       19465, 19483, 46977,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83702, 0, 3,
                                                                       83027, 46557, 83042,
                                                                       19483, 19501, 47007,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83747, 0, 3,
                                                                       83042, 46567, 83057,
                                                                       19501, 19519, 47037,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 83792, 0, 3,
                                                                       83072, 46587, 83117,
                                                                       19555, 19591, 47067,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 83882, 0, 3,
                                                                       83117, 46617, 83162,
                                                                       19591, 19627, 47127,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 83972, 0, 3,
                                                                       83162, 46647, 83207,
                                                                       19627, 19663, 47187,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84062, 0, 3,
                                                                       83207, 46677, 83252,
                                                                       19663, 19699, 47247,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84152, 0, 3,
                                                                       83252, 46707, 83297,
                                                                       19699, 19735, 47307,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84242, 0, 3,
                                                                       83297, 46737, 83342,
                                                                       19735, 19771, 47367,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84332, 0, 3,
                                                                       83342, 46767, 83387,
                                                                       19771, 19807, 47427,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84422, 0, 3,
                                                                       83387, 46797, 83432,
                                                                       19807, 19843, 47487,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84512, 0, 3,
                                                                       83432, 46827, 83477,
                                                                       19843, 19879, 47547,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84602, 0, 3,
                                                                       83477, 46857, 83522,
                                                                       19879, 19915, 47607,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84692, 0, 3,
                                                                       83522, 46887, 83567,
                                                                       19915, 19951, 47667,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84782, 0, 3,
                                                                       83567, 46917, 83612,
                                                                       19951, 19987, 47727,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84872, 0, 3,
                                                                       83612, 46947, 83657,
                                                                       19987, 20023, 47787,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84962, 0, 3,
                                                                       83657, 46977, 83702,
                                                                       20023, 20059, 47847,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 85052, 0, 3,
                                                                       83702, 47007, 83747,
                                                                       20059, 20095, 47907,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85142, 0, 3,
                                                                       83792, 47067, 83882,
                                                                       20167, 20227, 47967,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85292, 0, 3,
                                                                       83882, 47127, 83972,
                                                                       20227, 20287, 48067,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85442, 0, 3,
                                                                       83972, 47187, 84062,
                                                                       20287, 20347, 48167,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85592, 0, 3,
                                                                       84062, 47247, 84152,
                                                                       20347, 20407, 48267,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85742, 0, 3,
                                                                       84152, 47307, 84242,
                                                                       20407, 20467, 48367,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85892, 0, 3,
                                                                       84242, 47367, 84332,
                                                                       20467, 20527, 48467,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86042, 0, 3,
                                                                       84332, 47427, 84422,
                                                                       20527, 20587, 48567,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86192, 0, 3,
                                                                       84422, 47487, 84512,
                                                                       20587, 20647, 48667,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86342, 0, 3,
                                                                       84512, 47547, 84602,
                                                                       20647, 20707, 48767,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86492, 0, 3,
                                                                       84602, 47607, 84692,
                                                                       20707, 20767, 48867,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86642, 0, 3,
                                                                       84692, 47667, 84782,
                                                                       20767, 20827, 48967,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86792, 0, 3,
                                                                       84782, 47727, 84872,
                                                                       20827, 20887, 49067,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86942, 0, 3,
                                                                       84872, 47787, 84962,
                                                                       20887, 20947, 49167,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 87092, 0, 3,
                                                                       84962, 47847, 85052,
                                                                       20947, 21007, 49267,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 87242, 0, 3,
                                                                       85142, 47967, 85292,
                                                                       21127, 21217, 49367,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 87467, 0, 3,
                                                                       85292, 48067, 85442,
                                                                       21217, 21307, 49517,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 87692, 0, 3,
                                                                       85442, 48167, 85592,
                                                                       21307, 21397, 49667,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 87917, 0, 3,
                                                                       85592, 48267, 85742,
                                                                       21397, 21487, 49817,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 88142, 0, 3,
                                                                       85742, 48367, 85892,
                                                                       21487, 21577, 49967,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 88367, 0, 3,
                                                                       85892, 48467, 86042,
                                                                       21577, 21667, 50117,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 88592, 0, 3,
                                                                       86042, 48567, 86192,
                                                                       21667, 21757, 50267,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 88817, 0, 3,
                                                                       86192, 48667, 86342,
                                                                       21757, 21847, 50417,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 89042, 0, 3,
                                                                       86342, 48767, 86492,
                                                                       21847, 21937, 50567,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 89267, 0, 3,
                                                                       86492, 48867, 86642,
                                                                       21937, 22027, 50717,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 89492, 0, 3,
                                                                       86642, 48967, 86792,
                                                                       22027, 22117, 50867,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 89717, 0, 3,
                                                                       86792, 49067, 86942,
                                                                       22117, 22207, 51017,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 89942, 0, 3,
                                                                       86942, 49167, 87092,
                                                                       22207, 22297, 51167,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 90167, 0, 3,
                                                                       87242, 49367, 87467,
                                                                       22477, 22603, 51317,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 90482, 0, 3,
                                                                       87467, 49517, 87692,
                                                                       22603, 22729, 51527,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 90797, 0, 3,
                                                                       87692, 49667, 87917,
                                                                       22729, 22855, 51737,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 91112, 0, 3,
                                                                       87917, 49817, 88142,
                                                                       22855, 22981, 51947,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 91427, 0, 3,
                                                                       88142, 49967, 88367,
                                                                       22981, 23107, 52157,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 91742, 0, 3,
                                                                       88367, 50117, 88592,
                                                                       23107, 23233, 52367,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 92057, 0, 3,
                                                                       88592, 50267, 88817,
                                                                       23233, 23359, 52577,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 92372, 0, 3,
                                                                       88817, 50417, 89042,
                                                                       23359, 23485, 52787,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 92687, 0, 3,
                                                                       89042, 50567, 89267,
                                                                       23485, 23611, 52997,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 93002, 0, 3,
                                                                       89267, 50717, 89492,
                                                                       23611, 23737, 53207,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 93317, 0, 3,
                                                                       89492, 50867, 89717,
                                                                       23737, 23863, 53417,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 93632, 0, 3,
                                                                       89717, 51017, 89942,
                                                                       23863, 23989, 53627,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 93947, 0, 3,
                                                                       90167, 51317, 90482,
                                                                       24241, 24409, 53837,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 94367, 0, 3,
                                                                       90482, 51527, 90797,
                                                                       24409, 24577, 54117,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 94787, 0, 3,
                                                                       90797, 51737, 91112,
                                                                       24577, 24745, 54397,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 95207, 0, 3,
                                                                       91112, 51947, 91427,
                                                                       24745, 24913, 54677,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 95627, 0, 3,
                                                                       91427, 52157, 91742,
                                                                       24913, 25081, 54957,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 96047, 0, 3,
                                                                       91742, 52367, 92057,
                                                                       25081, 25249, 55237,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 96467, 0, 3,
                                                                       92057, 52577, 92372,
                                                                       25249, 25417, 55517,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 96887, 0, 3,
                                                                       92372, 52787, 92687,
                                                                       25417, 25585, 55797,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 97307, 0, 3,
                                                                       92687, 52997, 93002,
                                                                       25585, 25753, 56077,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 97727, 0, 3,
                                                                       93002, 53207, 93317,
                                                                       25753, 25921, 56357,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 98147, 0, 3,
                                                                       93317, 53417, 93632,
                                                                       25921, 26089, 56637,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 98567, 0, 3,
                                                                       93947, 53837, 94367,
                                                                       26425, 26641, 56917,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 99107, 0, 3,
                                                                       94367, 54117, 94787,
                                                                       26641, 26857, 57277,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 99647, 0, 3,
                                                                       94787, 54397, 95207,
                                                                       26857, 27073, 57637,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 100187, 0, 3,
                                                                       95207, 54677, 95627,
                                                                       27073, 27289, 57997,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 100727, 0, 3,
                                                                       95627, 54957, 96047,
                                                                       27289, 27505, 58357,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 101267, 0, 3,
                                                                       96047, 55237, 96467,
                                                                       27505, 27721, 58717,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 101807, 0, 3,
                                                                       96467, 55517, 96887,
                                                                       27721, 27937, 59077,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 102347, 0, 3,
                                                                       96887, 55797, 97307,
                                                                       27937, 28153, 59437,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 102887, 0, 3,
                                                                       97307, 56077, 97727,
                                                                       28153, 28369, 59797,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 103427, 0, 3,
                                                                       97727, 56357, 98147,
                                                                       28369, 28585, 60157,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 103967, 0, 3,
                                                                       98567, 56917, 99107,
                                                                       29017, 29287, 60517,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 104642, 0, 3,
                                                                       99107, 57277, 99647,
                                                                       29287, 29557, 60967,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 105317, 0, 3,
                                                                       99647, 57637, 100187,
                                                                       29557, 29827, 61417,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 105992, 0, 3,
                                                                       100187, 57997, 100727,
                                                                       29827, 30097, 61867,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 106667, 0, 3,
                                                                       100727, 58357, 101267,
                                                                       30097, 30367, 62317,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 107342, 0, 3,
                                                                       101267, 58717, 101807,
                                                                       30367, 30637, 62767,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 108017, 0, 3,
                                                                       101807, 59077, 102347,
                                                                       30637, 30907, 63217,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 108692, 0, 3,
                                                                       102347, 59437, 102887,
                                                                       30907, 31177, 63667,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 109367, 0, 3,
                                                                       102887, 59797, 103427,
                                                                       31177, 31447, 64117,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 110042, 0, 3,
                                                                       103967, 60517, 104642,
                                                                       31987, 32317, 64567,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 110867, 0, 3,
                                                                       104642, 60967, 105317,
                                                                       32317, 32647, 65117,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 111692, 0, 3,
                                                                       105317, 61417, 105992,
                                                                       32647, 32977, 65667,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 112517, 0, 3,
                                                                       105992, 61867, 106667,
                                                                       32977, 33307, 66217,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 113342, 0, 3,
                                                                       106667, 62317, 107342,
                                                                       33307, 33637, 66767,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 114167, 0, 3,
                                                                       107342, 62767, 108017,
                                                                       33637, 33967, 67317,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 114992, 0, 3,
                                                                       108017, 63217, 108692,
                                                                       33967, 34297, 67867,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 115817, 0, 3,
                                                                       108692, 63667, 109367,
                                                                       34297, 34627, 68417,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 116642, 0, 3,
                                                                       110042, 64567, 110867,
                                                                       35287, 35683, 68967,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 117632, 0, 3,
                                                                       110867, 65117, 111692,
                                                                       35683, 36079, 69627,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 118622, 0, 3,
                                                                       111692, 65667, 112517,
                                                                       36079, 36475, 70287,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 119612, 0, 3,
                                                                       112517, 66217, 113342,
                                                                       36475, 36871, 70947,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 120602, 0, 3,
                                                                       113342, 66767, 114167,
                                                                       36871, 37267, 71607,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 121592, 0, 3,
                                                                       114167, 67317, 114992,
                                                                       37267, 37663, 72267,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 122582, 0, 3,
                                                                       114992, 67867, 115817,
                                                                       37663, 38059, 72927,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 123572, 0, 3,
                                                                       116642, 68967, 117632,
                                                                       38851, 39319, 73587,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 124742, 0, 3,
                                                                       117632, 69627, 118622,
                                                                       39319, 39787, 74367,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 125912, 0, 3,
                                                                       118622, 70287, 119612,
                                                                       39787, 40255, 75147,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 127082, 0, 3,
                                                                       119612, 70947, 120602,
                                                                       40255, 40723, 75927,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 128252, 0, 3,
                                                                       120602, 71607, 121592,
                                                                       40723, 41191, 76707,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 129422, 0, 3,
                                                                       121592, 72267, 122582,
                                                                       41191, 41659, 77487,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 130592, 0, 3,
                                                                       123572, 73587, 124742,
                                                                       42595, 43141, 78267,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 131957, 0, 3,
                                                                       124742, 74367, 125912,
                                                                       43141, 43687, 79177,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 133322, 0, 3,
                                                                       125912, 75147, 127082,
                                                                       43687, 44233, 80087,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 134687, 0, 3,
                                                                       127082, 75927, 128252,
                                                                       44233, 44779, 80997,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 136052, 0, 3,
                                                                       128252, 76707, 129422,
                                                                       44779, 45325, 81907,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137417, 3, 46417,
                                                                       46427, 82847, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137438, 3, 46427,
                                                                       46437, 82862, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137459, 3, 46437,
                                                                       46447, 82877, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137480, 3, 46447,
                                                                       46457, 82892, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137501, 3, 46457,
                                                                       46467, 82907, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137522, 3, 46467,
                                                                       46477, 82922, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137543, 3, 46477,
                                                                       46487, 82937, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137564, 3, 46487,
                                                                       46497, 82952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137585, 3, 46497,
                                                                       46507, 82967, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137606, 3, 46507,
                                                                       46517, 82982, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137627, 3, 46517,
                                                                       46527, 82997, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137648, 3, 46527,
                                                                       46537, 83012, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137669, 3, 46537,
                                                                       46547, 83027, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137690, 3, 46547,
                                                                       46557, 83042, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137711, 3, 46557,
                                                                       46567, 83057, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 137732, 0, 3,
                                                                       137417, 82847, 137438,
                                                                       46587, 46617, 83162,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 137795, 0, 3,
                                                                       137438, 82862, 137459,
                                                                       46617, 46647, 83207,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 137858, 0, 3,
                                                                       137459, 82877, 137480,
                                                                       46647, 46677, 83252,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 137921, 0, 3,
                                                                       137480, 82892, 137501,
                                                                       46677, 46707, 83297,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 137984, 0, 3,
                                                                       137501, 82907, 137522,
                                                                       46707, 46737, 83342,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138047, 0, 3,
                                                                       137522, 82922, 137543,
                                                                       46737, 46767, 83387,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138110, 0, 3,
                                                                       137543, 82937, 137564,
                                                                       46767, 46797, 83432,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138173, 0, 3,
                                                                       137564, 82952, 137585,
                                                                       46797, 46827, 83477,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138236, 0, 3,
                                                                       137585, 82967, 137606,
                                                                       46827, 46857, 83522,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138299, 0, 3,
                                                                       137606, 82982, 137627,
                                                                       46857, 46887, 83567,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138362, 0, 3,
                                                                       137627, 82997, 137648,
                                                                       46887, 46917, 83612,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138425, 0, 3,
                                                                       137648, 83012, 137669,
                                                                       46917, 46947, 83657,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138488, 0, 3,
                                                                       137669, 83027, 137690,
                                                                       46947, 46977, 83702,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138551, 0, 3,
                                                                       137690, 83042, 137711,
                                                                       46977, 47007, 83747,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 138614, 0, 3,
                                                                       137732, 83162, 137795,
                                                                       47067, 47127, 83972,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 138740, 0, 3,
                                                                       137795, 83207, 137858,
                                                                       47127, 47187, 84062,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 138866, 0, 3,
                                                                       137858, 83252, 137921,
                                                                       47187, 47247, 84152,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 138992, 0, 3,
                                                                       137921, 83297, 137984,
                                                                       47247, 47307, 84242,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139118, 0, 3,
                                                                       137984, 83342, 138047,
                                                                       47307, 47367, 84332,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139244, 0, 3,
                                                                       138047, 83387, 138110,
                                                                       47367, 47427, 84422,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139370, 0, 3,
                                                                       138110, 83432, 138173,
                                                                       47427, 47487, 84512,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139496, 0, 3,
                                                                       138173, 83477, 138236,
                                                                       47487, 47547, 84602,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139622, 0, 3,
                                                                       138236, 83522, 138299,
                                                                       47547, 47607, 84692,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139748, 0, 3,
                                                                       138299, 83567, 138362,
                                                                       47607, 47667, 84782,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139874, 0, 3,
                                                                       138362, 83612, 138425,
                                                                       47667, 47727, 84872,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140000, 0, 3,
                                                                       138425, 83657, 138488,
                                                                       47727, 47787, 84962,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140126, 0, 3,
                                                                       138488, 83702, 138551,
                                                                       47787, 47847, 85052,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 140252, 0, 3,
                                                                       138614, 83972, 138740,
                                                                       47967, 48067, 85442,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 140462, 0, 3,
                                                                       138740, 84062, 138866,
                                                                       48067, 48167, 85592,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 140672, 0, 3,
                                                                       138866, 84152, 138992,
                                                                       48167, 48267, 85742,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 140882, 0, 3,
                                                                       138992, 84242, 139118,
                                                                       48267, 48367, 85892,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 141092, 0, 3,
                                                                       139118, 84332, 139244,
                                                                       48367, 48467, 86042,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 141302, 0, 3,
                                                                       139244, 84422, 139370,
                                                                       48467, 48567, 86192,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 141512, 0, 3,
                                                                       139370, 84512, 139496,
                                                                       48567, 48667, 86342,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 141722, 0, 3,
                                                                       139496, 84602, 139622,
                                                                       48667, 48767, 86492,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 141932, 0, 3,
                                                                       139622, 84692, 139748,
                                                                       48767, 48867, 86642,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 142142, 0, 3,
                                                                       139748, 84782, 139874,
                                                                       48867, 48967, 86792,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 142352, 0, 3,
                                                                       139874, 84872, 140000,
                                                                       48967, 49067, 86942,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 142562, 0, 3,
                                                                       140000, 84962, 140126,
                                                                       49067, 49167, 87092,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 142772, 0, 3,
                                                                       140252, 85442, 140462,
                                                                       49367, 49517, 87692,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 143087, 0, 3,
                                                                       140462, 85592, 140672,
                                                                       49517, 49667, 87917,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 143402, 0, 3,
                                                                       140672, 85742, 140882,
                                                                       49667, 49817, 88142,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 143717, 0, 3,
                                                                       140882, 85892, 141092,
                                                                       49817, 49967, 88367,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 144032, 0, 3,
                                                                       141092, 86042, 141302,
                                                                       49967, 50117, 88592,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 144347, 0, 3,
                                                                       141302, 86192, 141512,
                                                                       50117, 50267, 88817,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 144662, 0, 3,
                                                                       141512, 86342, 141722,
                                                                       50267, 50417, 89042,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 144977, 0, 3,
                                                                       141722, 86492, 141932,
                                                                       50417, 50567, 89267,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 145292, 0, 3,
                                                                       141932, 86642, 142142,
                                                                       50567, 50717, 89492,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 145607, 0, 3,
                                                                       142142, 86792, 142352,
                                                                       50717, 50867, 89717,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 145922, 0, 3,
                                                                       142352, 86942, 142562,
                                                                       50867, 51017, 89942,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 146237, 0, 3,
                                                                       142772, 87692, 143087,
                                                                       51317, 51527, 90797,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 146678, 0, 3,
                                                                       143087, 87917, 143402,
                                                                       51527, 51737, 91112,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 147119, 0, 3,
                                                                       143402, 88142, 143717,
                                                                       51737, 51947, 91427,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 147560, 0, 3,
                                                                       143717, 88367, 144032,
                                                                       51947, 52157, 91742,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 148001, 0, 3,
                                                                       144032, 88592, 144347,
                                                                       52157, 52367, 92057,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 148442, 0, 3,
                                                                       144347, 88817, 144662,
                                                                       52367, 52577, 92372,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 148883, 0, 3,
                                                                       144662, 89042, 144977,
                                                                       52577, 52787, 92687,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 149324, 0, 3,
                                                                       144977, 89267, 145292,
                                                                       52787, 52997, 93002,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 149765, 0, 3,
                                                                       145292, 89492, 145607,
                                                                       52997, 53207, 93317,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 150206, 0, 3,
                                                                       145607, 89717, 145922,
                                                                       53207, 53417, 93632,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 150647, 0, 3,
                                                                       146237, 90797, 146678,
                                                                       53837, 54117, 94787,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 151235, 0, 3,
                                                                       146678, 91112, 147119,
                                                                       54117, 54397, 95207,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 151823, 0, 3,
                                                                       147119, 91427, 147560,
                                                                       54397, 54677, 95627,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 152411, 0, 3,
                                                                       147560, 91742, 148001,
                                                                       54677, 54957, 96047,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 152999, 0, 3,
                                                                       148001, 92057, 148442,
                                                                       54957, 55237, 96467,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 153587, 0, 3,
                                                                       148442, 92372, 148883,
                                                                       55237, 55517, 96887,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 154175, 0, 3,
                                                                       148883, 92687, 149324,
                                                                       55517, 55797, 97307,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 154763, 0, 3,
                                                                       149324, 93002, 149765,
                                                                       55797, 56077, 97727,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 155351, 0, 3,
                                                                       149765, 93317, 150206,
                                                                       56077, 56357, 98147,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 155939, 0, 3,
                                                                       150647, 94787, 151235,
                                                                       56917, 57277, 99647,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 156695, 0, 3,
                                                                       151235, 95207, 151823,
                                                                       57277, 57637, 100187,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 157451, 0, 3,
                                                                       151823, 95627, 152411,
                                                                       57637, 57997, 100727,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 158207, 0, 3,
                                                                       152411, 96047, 152999,
                                                                       57997, 58357, 101267,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 158963, 0, 3,
                                                                       152999, 96467, 153587,
                                                                       58357, 58717, 101807,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 159719, 0, 3,
                                                                       153587, 96887, 154175,
                                                                       58717, 59077, 102347,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 160475, 0, 3,
                                                                       154175, 97307, 154763,
                                                                       59077, 59437, 102887,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 161231, 0, 3,
                                                                       154763, 97727, 155351,
                                                                       59437, 59797, 103427,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 161987, 0, 3,
                                                                       155939, 99647, 156695,
                                                                       60517, 60967, 105317,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 162932, 0, 3,
                                                                       156695, 100187, 157451,
                                                                       60967, 61417, 105992,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 163877, 0, 3,
                                                                       157451, 100727, 158207,
                                                                       61417, 61867, 106667,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 164822, 0, 3,
                                                                       158207, 101267, 158963,
                                                                       61867, 62317, 107342,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 165767, 0, 3,
                                                                       158963, 101807, 159719,
                                                                       62317, 62767, 108017,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 166712, 0, 3,
                                                                       159719, 102347, 160475,
                                                                       62767, 63217, 108692,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 167657, 0, 3,
                                                                       160475, 102887, 161231,
                                                                       63217, 63667, 109367,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 168602, 0, 3,
                                                                       161987, 105317, 162932,
                                                                       64567, 65117, 111692,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 169757, 0, 3,
                                                                       162932, 105992, 163877,
                                                                       65117, 65667, 112517,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 170912, 0, 3,
                                                                       163877, 106667, 164822,
                                                                       65667, 66217, 113342,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 172067, 0, 3,
                                                                       164822, 107342, 165767,
                                                                       66217, 66767, 114167,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 173222, 0, 3,
                                                                       165767, 108017, 166712,
                                                                       66767, 67317, 114992,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 174377, 0, 3,
                                                                       166712, 108692, 167657,
                                                                       67317, 67867, 115817,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 175532, 0, 3,
                                                                       168602, 111692, 169757,
                                                                       68967, 69627, 118622,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 176918, 0, 3,
                                                                       169757, 112517, 170912,
                                                                       69627, 70287, 119612,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 178304, 0, 3,
                                                                       170912, 113342, 172067,
                                                                       70287, 70947, 120602,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 179690, 0, 3,
                                                                       172067, 114167, 173222,
                                                                       70947, 71607, 121592,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 181076, 0, 3,
                                                                       173222, 114992, 174377,
                                                                       71607, 72267, 122582,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 182462, 0, 3,
                                                                       175532, 118622, 176918,
                                                                       73587, 74367, 125912,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 184100, 0, 3,
                                                                       176918, 119612, 178304,
                                                                       74367, 75147, 127082,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 185738, 0, 3,
                                                                       178304, 120602, 179690,
                                                                       75147, 75927, 128252,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 187376, 0, 3,
                                                                       179690, 121592, 181076,
                                                                       75927, 76707, 129422,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 189014, 0, 3,
                                                                       182462, 125912, 184100,
                                                                       78267, 79177, 133322,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 190925, 0, 3,
                                                                       184100, 127082, 185738,
                                                                       79177, 80087, 134687,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 192836, 0, 3,
                                                                       185738, 128252, 187376,
                                                                       80087, 80997, 136052,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194747, 3, 82817,
                                                                       82832, 137417, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194775, 3, 82832,
                                                                       82847, 137438, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194803, 3, 82847,
                                                                       82862, 137459, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194831, 3, 82862,
                                                                       82877, 137480, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194859, 3, 82877,
                                                                       82892, 137501, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194887, 3, 82892,
                                                                       82907, 137522, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194915, 3, 82907,
                                                                       82922, 137543, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194943, 3, 82922,
                                                                       82937, 137564, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194971, 3, 82937,
                                                                       82952, 137585, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194999, 3, 82952,
                                                                       82967, 137606, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 195027, 3, 82967,
                                                                       82982, 137627, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 195055, 3, 82982,
                                                                       82997, 137648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 195083, 3, 82997,
                                                                       83012, 137669, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 195111, 3, 83012,
                                                                       83027, 137690, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 195139, 3, 83027,
                                                                       83042, 137711, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195167, 0, 3,
                                                                       194747, 137417, 194775,
                                                                       83072, 83117, 137732,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195251, 0, 3,
                                                                       194775, 137438, 194803,
                                                                       83117, 83162, 137795,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195335, 0, 3,
                                                                       194803, 137459, 194831,
                                                                       83162, 83207, 137858,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195419, 0, 3,
                                                                       194831, 137480, 194859,
                                                                       83207, 83252, 137921,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195503, 0, 3,
                                                                       194859, 137501, 194887,
                                                                       83252, 83297, 137984,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195587, 0, 3,
                                                                       194887, 137522, 194915,
                                                                       83297, 83342, 138047,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195671, 0, 3,
                                                                       194915, 137543, 194943,
                                                                       83342, 83387, 138110,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195755, 0, 3,
                                                                       194943, 137564, 194971,
                                                                       83387, 83432, 138173,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195839, 0, 3,
                                                                       194971, 137585, 194999,
                                                                       83432, 83477, 138236,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195923, 0, 3,
                                                                       194999, 137606, 195027,
                                                                       83477, 83522, 138299,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 196007, 0, 3,
                                                                       195027, 137627, 195055,
                                                                       83522, 83567, 138362,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 196091, 0, 3,
                                                                       195055, 137648, 195083,
                                                                       83567, 83612, 138425,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 196175, 0, 3,
                                                                       195083, 137669, 195111,
                                                                       83612, 83657, 138488,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 196259, 0, 3,
                                                                       195111, 137690, 195139,
                                                                       83657, 83702, 138551,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 196343, 0, 3,
                                                                       195167, 137732, 195251,
                                                                       83792, 83882, 138614,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 196511, 0, 3,
                                                                       195251, 137795, 195335,
                                                                       83882, 83972, 138740,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 196679, 0, 3,
                                                                       195335, 137858, 195419,
                                                                       83972, 84062, 138866,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 196847, 0, 3,
                                                                       195419, 137921, 195503,
                                                                       84062, 84152, 138992,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197015, 0, 3,
                                                                       195503, 137984, 195587,
                                                                       84152, 84242, 139118,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197183, 0, 3,
                                                                       195587, 138047, 195671,
                                                                       84242, 84332, 139244,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197351, 0, 3,
                                                                       195671, 138110, 195755,
                                                                       84332, 84422, 139370,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197519, 0, 3,
                                                                       195755, 138173, 195839,
                                                                       84422, 84512, 139496,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197687, 0, 3,
                                                                       195839, 138236, 195923,
                                                                       84512, 84602, 139622,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197855, 0, 3,
                                                                       195923, 138299, 196007,
                                                                       84602, 84692, 139748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 198023, 0, 3,
                                                                       196007, 138362, 196091,
                                                                       84692, 84782, 139874,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 198191, 0, 3,
                                                                       196091, 138425, 196175,
                                                                       84782, 84872, 140000,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 198359, 0, 3,
                                                                       196175, 138488, 196259,
                                                                       84872, 84962, 140126,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 198527, 0, 3,
                                                                       196343, 138614, 196511,
                                                                       85142, 85292, 140252,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 198807, 0, 3,
                                                                       196511, 138740, 196679,
                                                                       85292, 85442, 140462,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 199087, 0, 3,
                                                                       196679, 138866, 196847,
                                                                       85442, 85592, 140672,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 199367, 0, 3,
                                                                       196847, 138992, 197015,
                                                                       85592, 85742, 140882,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 199647, 0, 3,
                                                                       197015, 139118, 197183,
                                                                       85742, 85892, 141092,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 199927, 0, 3,
                                                                       197183, 139244, 197351,
                                                                       85892, 86042, 141302,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 200207, 0, 3,
                                                                       197351, 139370, 197519,
                                                                       86042, 86192, 141512,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 200487, 0, 3,
                                                                       197519, 139496, 197687,
                                                                       86192, 86342, 141722,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 200767, 0, 3,
                                                                       197687, 139622, 197855,
                                                                       86342, 86492, 141932,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 201047, 0, 3,
                                                                       197855, 139748, 198023,
                                                                       86492, 86642, 142142,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 201327, 0, 3,
                                                                       198023, 139874, 198191,
                                                                       86642, 86792, 142352,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 201607, 0, 3,
                                                                       198191, 140000, 198359,
                                                                       86792, 86942, 142562,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 201887, 0, 3,
                                                                       198527, 140252, 198807,
                                                                       87242, 87467, 142772,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 202307, 0, 3,
                                                                       198807, 140462, 199087,
                                                                       87467, 87692, 143087,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 202727, 0, 3,
                                                                       199087, 140672, 199367,
                                                                       87692, 87917, 143402,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 203147, 0, 3,
                                                                       199367, 140882, 199647,
                                                                       87917, 88142, 143717,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 203567, 0, 3,
                                                                       199647, 141092, 199927,
                                                                       88142, 88367, 144032,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 203987, 0, 3,
                                                                       199927, 141302, 200207,
                                                                       88367, 88592, 144347,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 204407, 0, 3,
                                                                       200207, 141512, 200487,
                                                                       88592, 88817, 144662,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 204827, 0, 3,
                                                                       200487, 141722, 200767,
                                                                       88817, 89042, 144977,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 205247, 0, 3,
                                                                       200767, 141932, 201047,
                                                                       89042, 89267, 145292,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 205667, 0, 3,
                                                                       201047, 142142, 201327,
                                                                       89267, 89492, 145607,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 206087, 0, 3,
                                                                       201327, 142352, 201607,
                                                                       89492, 89717, 145922,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 206507, 0, 3,
                                                                       201887, 142772, 202307,
                                                                       90167, 90482, 146237,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 207095, 0, 3,
                                                                       202307, 143087, 202727,
                                                                       90482, 90797, 146678,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 207683, 0, 3,
                                                                       202727, 143402, 203147,
                                                                       90797, 91112, 147119,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 208271, 0, 3,
                                                                       203147, 143717, 203567,
                                                                       91112, 91427, 147560,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 208859, 0, 3,
                                                                       203567, 144032, 203987,
                                                                       91427, 91742, 148001,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 209447, 0, 3,
                                                                       203987, 144347, 204407,
                                                                       91742, 92057, 148442,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 210035, 0, 3,
                                                                       204407, 144662, 204827,
                                                                       92057, 92372, 148883,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 210623, 0, 3,
                                                                       204827, 144977, 205247,
                                                                       92372, 92687, 149324,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 211211, 0, 3,
                                                                       205247, 145292, 205667,
                                                                       92687, 93002, 149765,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 211799, 0, 3,
                                                                       205667, 145607, 206087,
                                                                       93002, 93317, 150206,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 212387, 0, 3,
                                                                       206507, 146237, 207095,
                                                                       93947, 94367, 150647,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 213171, 0, 3,
                                                                       207095, 146678, 207683,
                                                                       94367, 94787, 151235,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 213955, 0, 3,
                                                                       207683, 147119, 208271,
                                                                       94787, 95207, 151823,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 214739, 0, 3,
                                                                       208271, 147560, 208859,
                                                                       95207, 95627, 152411,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 215523, 0, 3,
                                                                       208859, 148001, 209447,
                                                                       95627, 96047, 152999,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 216307, 0, 3,
                                                                       209447, 148442, 210035,
                                                                       96047, 96467, 153587,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 217091, 0, 3,
                                                                       210035, 148883, 210623,
                                                                       96467, 96887, 154175,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 217875, 0, 3,
                                                                       210623, 149324, 211211,
                                                                       96887, 97307, 154763,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 218659, 0, 3,
                                                                       211211, 149765, 211799,
                                                                       97307, 97727, 155351,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 219443, 0, 3,
                                                                       212387, 150647, 213171,
                                                                       98567, 99107, 155939,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 220451, 0, 3,
                                                                       213171, 151235, 213955,
                                                                       99107, 99647, 156695,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 221459, 0, 3,
                                                                       213955, 151823, 214739,
                                                                       99647, 100187, 157451,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 222467, 0, 3,
                                                                       214739, 152411, 215523,
                                                                       100187, 100727, 158207,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 223475, 0, 3,
                                                                       215523, 152999, 216307,
                                                                       100727, 101267, 158963,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 224483, 0, 3,
                                                                       216307, 153587, 217091,
                                                                       101267, 101807, 159719,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 225491, 0, 3,
                                                                       217091, 154175, 217875,
                                                                       101807, 102347, 160475,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 226499, 0, 3,
                                                                       217875, 154763, 218659,
                                                                       102347, 102887, 161231,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 227507, 0, 3,
                                                                       219443, 155939, 220451,
                                                                       103967, 104642, 161987,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 228767, 0, 3,
                                                                       220451, 156695, 221459,
                                                                       104642, 105317, 162932,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 230027, 0, 3,
                                                                       221459, 157451, 222467,
                                                                       105317, 105992, 163877,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 231287, 0, 3,
                                                                       222467, 158207, 223475,
                                                                       105992, 106667, 164822,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 232547, 0, 3,
                                                                       223475, 158963, 224483,
                                                                       106667, 107342, 165767,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 233807, 0, 3,
                                                                       224483, 159719, 225491,
                                                                       107342, 108017, 166712,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 235067, 0, 3,
                                                                       225491, 160475, 226499,
                                                                       108017, 108692, 167657,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 236327, 0, 3,
                                                                       227507, 161987, 228767,
                                                                       110042, 110867, 168602,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 237867, 0, 3,
                                                                       228767, 162932, 230027,
                                                                       110867, 111692, 169757,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 239407, 0, 3,
                                                                       230027, 163877, 231287,
                                                                       111692, 112517, 170912,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 240947, 0, 3,
                                                                       231287, 164822, 232547,
                                                                       112517, 113342, 172067,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 242487, 0, 3,
                                                                       232547, 165767, 233807,
                                                                       113342, 114167, 173222,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 244027, 0, 3,
                                                                       233807, 166712, 235067,
                                                                       114167, 114992, 174377,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 245567, 0, 3,
                                                                       236327, 168602, 237867,
                                                                       116642, 117632, 175532,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 247415, 0, 3,
                                                                       237867, 169757, 239407,
                                                                       117632, 118622, 176918,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 249263, 0, 3,
                                                                       239407, 170912, 240947,
                                                                       118622, 119612, 178304,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 251111, 0, 3,
                                                                       240947, 172067, 242487,
                                                                       119612, 120602, 179690,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 252959, 0, 3,
                                                                       242487, 173222, 244027,
                                                                       120602, 121592, 181076,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 254807, 0, 3,
                                                                       245567, 175532, 247415,
                                                                       123572, 124742, 182462,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 256991, 0, 3,
                                                                       247415, 176918, 249263,
                                                                       124742, 125912, 184100,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 259175, 0, 3,
                                                                       249263, 178304, 251111,
                                                                       125912, 127082, 185738,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 261359, 0, 3,
                                                                       251111, 179690, 252959,
                                                                       127082, 128252, 187376,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 263543, 0, 3,
                                                                       254807, 182462, 256991,
                                                                       130592, 131957, 189014,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 266091, 0, 3,
                                                                       256991, 184100, 259175,
                                                                       131957, 133322, 190925,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 268639, 0, 3,
                                                                       259175, 185738, 261359,
                                                                       133322, 134687, 192836,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271187, 3, 137417,
                                                                       137438, 194803, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271223, 3, 137438,
                                                                       137459, 194831, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271259, 3, 137459,
                                                                       137480, 194859, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271295, 3, 137480,
                                                                       137501, 194887, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271331, 3, 137501,
                                                                       137522, 194915, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271367, 3, 137522,
                                                                       137543, 194943, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271403, 3, 137543,
                                                                       137564, 194971, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271439, 3, 137564,
                                                                       137585, 194999, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271475, 3, 137585,
                                                                       137606, 195027, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271511, 3, 137606,
                                                                       137627, 195055, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271547, 3, 137627,
                                                                       137648, 195083, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271583, 3, 137648,
                                                                       137669, 195111, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271619, 3, 137669,
                                                                       137690, 195139, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 271655, 0, 3,
                                                                       271187, 194803, 271223,
                                                                       137732, 137795, 195335,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 271763, 0, 3,
                                                                       271223, 194831, 271259,
                                                                       137795, 137858, 195419,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 271871, 0, 3,
                                                                       271259, 194859, 271295,
                                                                       137858, 137921, 195503,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 271979, 0, 3,
                                                                       271295, 194887, 271331,
                                                                       137921, 137984, 195587,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272087, 0, 3,
                                                                       271331, 194915, 271367,
                                                                       137984, 138047, 195671,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272195, 0, 3,
                                                                       271367, 194943, 271403,
                                                                       138047, 138110, 195755,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272303, 0, 3,
                                                                       271403, 194971, 271439,
                                                                       138110, 138173, 195839,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272411, 0, 3,
                                                                       271439, 194999, 271475,
                                                                       138173, 138236, 195923,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272519, 0, 3,
                                                                       271475, 195027, 271511,
                                                                       138236, 138299, 196007,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272627, 0, 3,
                                                                       271511, 195055, 271547,
                                                                       138299, 138362, 196091,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272735, 0, 3,
                                                                       271547, 195083, 271583,
                                                                       138362, 138425, 196175,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272843, 0, 3,
                                                                       271583, 195111, 271619,
                                                                       138425, 138488, 196259,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 272951, 0, 3,
                                                                       271655, 195335, 271763,
                                                                       138614, 138740, 196679,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 273167, 0, 3,
                                                                       271763, 195419, 271871,
                                                                       138740, 138866, 196847,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 273383, 0, 3,
                                                                       271871, 195503, 271979,
                                                                       138866, 138992, 197015,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 273599, 0, 3,
                                                                       271979, 195587, 272087,
                                                                       138992, 139118, 197183,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 273815, 0, 3,
                                                                       272087, 195671, 272195,
                                                                       139118, 139244, 197351,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 274031, 0, 3,
                                                                       272195, 195755, 272303,
                                                                       139244, 139370, 197519,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 274247, 0, 3,
                                                                       272303, 195839, 272411,
                                                                       139370, 139496, 197687,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 274463, 0, 3,
                                                                       272411, 195923, 272519,
                                                                       139496, 139622, 197855,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 274679, 0, 3,
                                                                       272519, 196007, 272627,
                                                                       139622, 139748, 198023,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 274895, 0, 3,
                                                                       272627, 196091, 272735,
                                                                       139748, 139874, 198191,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 275111, 0, 3,
                                                                       272735, 196175, 272843,
                                                                       139874, 140000, 198359,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 275327, 0, 3,
                                                                       272951, 196679, 273167,
                                                                       140252, 140462, 199087,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 275687, 0, 3,
                                                                       273167, 196847, 273383,
                                                                       140462, 140672, 199367,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 276047, 0, 3,
                                                                       273383, 197015, 273599,
                                                                       140672, 140882, 199647,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 276407, 0, 3,
                                                                       273599, 197183, 273815,
                                                                       140882, 141092, 199927,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 276767, 0, 3,
                                                                       273815, 197351, 274031,
                                                                       141092, 141302, 200207,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 277127, 0, 3,
                                                                       274031, 197519, 274247,
                                                                       141302, 141512, 200487,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 277487, 0, 3,
                                                                       274247, 197687, 274463,
                                                                       141512, 141722, 200767,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 277847, 0, 3,
                                                                       274463, 197855, 274679,
                                                                       141722, 141932, 201047,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 278207, 0, 3,
                                                                       274679, 198023, 274895,
                                                                       141932, 142142, 201327,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 278567, 0, 3,
                                                                       274895, 198191, 275111,
                                                                       142142, 142352, 201607,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 278927, 0, 3,
                                                                       275327, 199087, 275687,
                                                                       142772, 143087, 202727,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 279467, 0, 3,
                                                                       275687, 199367, 276047,
                                                                       143087, 143402, 203147,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 280007, 0, 3,
                                                                       276047, 199647, 276407,
                                                                       143402, 143717, 203567,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 280547, 0, 3,
                                                                       276407, 199927, 276767,
                                                                       143717, 144032, 203987,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 281087, 0, 3,
                                                                       276767, 200207, 277127,
                                                                       144032, 144347, 204407,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 281627, 0, 3,
                                                                       277127, 200487, 277487,
                                                                       144347, 144662, 204827,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 282167, 0, 3,
                                                                       277487, 200767, 277847,
                                                                       144662, 144977, 205247,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 282707, 0, 3,
                                                                       277847, 201047, 278207,
                                                                       144977, 145292, 205667,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 283247, 0, 3,
                                                                       278207, 201327, 278567,
                                                                       145292, 145607, 206087,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 283787, 0, 3,
                                                                       278927, 202727, 279467,
                                                                       146237, 146678, 207683,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 284543, 0, 3,
                                                                       279467, 203147, 280007,
                                                                       146678, 147119, 208271,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 285299, 0, 3,
                                                                       280007, 203567, 280547,
                                                                       147119, 147560, 208859,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 286055, 0, 3,
                                                                       280547, 203987, 281087,
                                                                       147560, 148001, 209447,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 286811, 0, 3,
                                                                       281087, 204407, 281627,
                                                                       148001, 148442, 210035,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 287567, 0, 3,
                                                                       281627, 204827, 282167,
                                                                       148442, 148883, 210623,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 288323, 0, 3,
                                                                       282167, 205247, 282707,
                                                                       148883, 149324, 211211,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 289079, 0, 3,
                                                                       282707, 205667, 283247,
                                                                       149324, 149765, 211799,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 289835, 0, 3,
                                                                       283787, 207683, 284543,
                                                                       150647, 151235, 213955,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 290843, 0, 3,
                                                                       284543, 208271, 285299,
                                                                       151235, 151823, 214739,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 291851, 0, 3,
                                                                       285299, 208859, 286055,
                                                                       151823, 152411, 215523,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 292859, 0, 3,
                                                                       286055, 209447, 286811,
                                                                       152411, 152999, 216307,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 293867, 0, 3,
                                                                       286811, 210035, 287567,
                                                                       152999, 153587, 217091,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 294875, 0, 3,
                                                                       287567, 210623, 288323,
                                                                       153587, 154175, 217875,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 295883, 0, 3,
                                                                       288323, 211211, 289079,
                                                                       154175, 154763, 218659,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 296891, 0, 3,
                                                                       289835, 213955, 290843,
                                                                       155939, 156695, 221459,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 298187, 0, 3,
                                                                       290843, 214739, 291851,
                                                                       156695, 157451, 222467,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 299483, 0, 3,
                                                                       291851, 215523, 292859,
                                                                       157451, 158207, 223475,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 300779, 0, 3,
                                                                       292859, 216307, 293867,
                                                                       158207, 158963, 224483,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 302075, 0, 3,
                                                                       293867, 217091, 294875,
                                                                       158963, 159719, 225491,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 303371, 0, 3,
                                                                       294875, 217875, 295883,
                                                                       159719, 160475, 226499,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 304667, 0, 3,
                                                                       296891, 221459, 298187,
                                                                       161987, 162932, 230027,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 306287, 0, 3,
                                                                       298187, 222467, 299483,
                                                                       162932, 163877, 231287,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 307907, 0, 3,
                                                                       299483, 223475, 300779,
                                                                       163877, 164822, 232547,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 309527, 0, 3,
                                                                       300779, 224483, 302075,
                                                                       164822, 165767, 233807,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 311147, 0, 3,
                                                                       302075, 225491, 303371,
                                                                       165767, 166712, 235067,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 312767, 0, 3,
                                                                       304667, 230027, 306287,
                                                                       168602, 169757, 239407,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 314747, 0, 3,
                                                                       306287, 231287, 307907,
                                                                       169757, 170912, 240947,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 316727, 0, 3,
                                                                       307907, 232547, 309527,
                                                                       170912, 172067, 242487,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 318707, 0, 3,
                                                                       309527, 233807, 311147,
                                                                       172067, 173222, 244027,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 320687, 0, 3,
                                                                       312767, 239407, 314747,
                                                                       175532, 176918, 249263,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 323063, 0, 3,
                                                                       314747, 240947, 316727,
                                                                       176918, 178304, 251111,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 325439, 0, 3,
                                                                       316727, 242487, 318707,
                                                                       178304, 179690, 252959,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 327815, 0, 3,
                                                                       320687, 249263, 323063,
                                                                       182462, 184100, 259175,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 330623, 0, 3,
                                                                       323063, 251111, 325439,
                                                                       184100, 185738, 261359,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsk_three_center_electron_repulsion_0(buffer, 333431, 0, 3,
                                                                       327815, 259175, 330623,
                                                                       189014, 190925, 268639,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336707, 3, 194747,
                                                                       194775, 271187, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336752, 3, 194775,
                                                                       194803, 271223, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336797, 3, 194803,
                                                                       194831, 271259, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336842, 3, 194831,
                                                                       194859, 271295, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336887, 3, 194859,
                                                                       194887, 271331, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336932, 3, 194887,
                                                                       194915, 271367, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336977, 3, 194915,
                                                                       194943, 271403, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337022, 3, 194943,
                                                                       194971, 271439, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337067, 3, 194971,
                                                                       194999, 271475, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337112, 3, 194999,
                                                                       195027, 271511, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337157, 3, 195027,
                                                                       195055, 271547, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337202, 3, 195055,
                                                                       195083, 271583, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337247, 3, 195083,
                                                                       195111, 271619, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337292, 0, 3,
                                                                       336707, 271187, 336752,
                                                                       195167, 195251, 271655,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337427, 0, 3,
                                                                       336752, 271223, 336797,
                                                                       195251, 195335, 271763,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337562, 0, 3,
                                                                       336797, 271259, 336842,
                                                                       195335, 195419, 271871,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337697, 0, 3,
                                                                       336842, 271295, 336887,
                                                                       195419, 195503, 271979,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337832, 0, 3,
                                                                       336887, 271331, 336932,
                                                                       195503, 195587, 272087,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337967, 0, 3,
                                                                       336932, 271367, 336977,
                                                                       195587, 195671, 272195,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338102, 0, 3,
                                                                       336977, 271403, 337022,
                                                                       195671, 195755, 272303,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338237, 0, 3,
                                                                       337022, 271439, 337067,
                                                                       195755, 195839, 272411,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338372, 0, 3,
                                                                       337067, 271475, 337112,
                                                                       195839, 195923, 272519,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338507, 0, 3,
                                                                       337112, 271511, 337157,
                                                                       195923, 196007, 272627,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338642, 0, 3,
                                                                       337157, 271547, 337202,
                                                                       196007, 196091, 272735,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338777, 0, 3,
                                                                       337202, 271583, 337247,
                                                                       196091, 196175, 272843,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 338912, 0, 3,
                                                                       337292, 271655, 337427,
                                                                       196343, 196511, 272951,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 339182, 0, 3,
                                                                       337427, 271763, 337562,
                                                                       196511, 196679, 273167,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 339452, 0, 3,
                                                                       337562, 271871, 337697,
                                                                       196679, 196847, 273383,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 339722, 0, 3,
                                                                       337697, 271979, 337832,
                                                                       196847, 197015, 273599,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 339992, 0, 3,
                                                                       337832, 272087, 337967,
                                                                       197015, 197183, 273815,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 340262, 0, 3,
                                                                       337967, 272195, 338102,
                                                                       197183, 197351, 274031,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 340532, 0, 3,
                                                                       338102, 272303, 338237,
                                                                       197351, 197519, 274247,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 340802, 0, 3,
                                                                       338237, 272411, 338372,
                                                                       197519, 197687, 274463,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 341072, 0, 3,
                                                                       338372, 272519, 338507,
                                                                       197687, 197855, 274679,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 341342, 0, 3,
                                                                       338507, 272627, 338642,
                                                                       197855, 198023, 274895,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 341612, 0, 3,
                                                                       338642, 272735, 338777,
                                                                       198023, 198191, 275111,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 341882, 0, 3,
                                                                       338912, 272951, 339182,
                                                                       198527, 198807, 275327,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 342332, 0, 3,
                                                                       339182, 273167, 339452,
                                                                       198807, 199087, 275687,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 342782, 0, 3,
                                                                       339452, 273383, 339722,
                                                                       199087, 199367, 276047,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 343232, 0, 3,
                                                                       339722, 273599, 339992,
                                                                       199367, 199647, 276407,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 343682, 0, 3,
                                                                       339992, 273815, 340262,
                                                                       199647, 199927, 276767,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 344132, 0, 3,
                                                                       340262, 274031, 340532,
                                                                       199927, 200207, 277127,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 344582, 0, 3,
                                                                       340532, 274247, 340802,
                                                                       200207, 200487, 277487,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 345032, 0, 3,
                                                                       340802, 274463, 341072,
                                                                       200487, 200767, 277847,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 345482, 0, 3,
                                                                       341072, 274679, 341342,
                                                                       200767, 201047, 278207,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 345932, 0, 3,
                                                                       341342, 274895, 341612,
                                                                       201047, 201327, 278567,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 346382, 0, 3,
                                                                       341882, 275327, 342332,
                                                                       201887, 202307, 278927,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 347057, 0, 3,
                                                                       342332, 275687, 342782,
                                                                       202307, 202727, 279467,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 347732, 0, 3,
                                                                       342782, 276047, 343232,
                                                                       202727, 203147, 280007,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 348407, 0, 3,
                                                                       343232, 276407, 343682,
                                                                       203147, 203567, 280547,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 349082, 0, 3,
                                                                       343682, 276767, 344132,
                                                                       203567, 203987, 281087,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 349757, 0, 3,
                                                                       344132, 277127, 344582,
                                                                       203987, 204407, 281627,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 350432, 0, 3,
                                                                       344582, 277487, 345032,
                                                                       204407, 204827, 282167,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 351107, 0, 3,
                                                                       345032, 277847, 345482,
                                                                       204827, 205247, 282707,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 351782, 0, 3,
                                                                       345482, 278207, 345932,
                                                                       205247, 205667, 283247,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 352457, 0, 3,
                                                                       346382, 278927, 347057,
                                                                       206507, 207095, 283787,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 353402, 0, 3,
                                                                       347057, 279467, 347732,
                                                                       207095, 207683, 284543,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 354347, 0, 3,
                                                                       347732, 280007, 348407,
                                                                       207683, 208271, 285299,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 355292, 0, 3,
                                                                       348407, 280547, 349082,
                                                                       208271, 208859, 286055,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 356237, 0, 3,
                                                                       349082, 281087, 349757,
                                                                       208859, 209447, 286811,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 357182, 0, 3,
                                                                       349757, 281627, 350432,
                                                                       209447, 210035, 287567,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 358127, 0, 3,
                                                                       350432, 282167, 351107,
                                                                       210035, 210623, 288323,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 359072, 0, 3,
                                                                       351107, 282707, 351782,
                                                                       210623, 211211, 289079,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 360017, 0, 3,
                                                                       352457, 283787, 353402,
                                                                       212387, 213171, 289835,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 361277, 0, 3,
                                                                       353402, 284543, 354347,
                                                                       213171, 213955, 290843,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 362537, 0, 3,
                                                                       354347, 285299, 355292,
                                                                       213955, 214739, 291851,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 363797, 0, 3,
                                                                       355292, 286055, 356237,
                                                                       214739, 215523, 292859,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 365057, 0, 3,
                                                                       356237, 286811, 357182,
                                                                       215523, 216307, 293867,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 366317, 0, 3,
                                                                       357182, 287567, 358127,
                                                                       216307, 217091, 294875,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 367577, 0, 3,
                                                                       358127, 288323, 359072,
                                                                       217091, 217875, 295883,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 368837, 0, 3,
                                                                       360017, 289835, 361277,
                                                                       219443, 220451, 296891,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 370457, 0, 3,
                                                                       361277, 290843, 362537,
                                                                       220451, 221459, 298187,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 372077, 0, 3,
                                                                       362537, 291851, 363797,
                                                                       221459, 222467, 299483,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 373697, 0, 3,
                                                                       363797, 292859, 365057,
                                                                       222467, 223475, 300779,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 375317, 0, 3,
                                                                       365057, 293867, 366317,
                                                                       223475, 224483, 302075,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 376937, 0, 3,
                                                                       366317, 294875, 367577,
                                                                       224483, 225491, 303371,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 378557, 0, 3,
                                                                       368837, 296891, 370457,
                                                                       227507, 228767, 304667,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 380582, 0, 3,
                                                                       370457, 298187, 372077,
                                                                       228767, 230027, 306287,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 382607, 0, 3,
                                                                       372077, 299483, 373697,
                                                                       230027, 231287, 307907,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 384632, 0, 3,
                                                                       373697, 300779, 375317,
                                                                       231287, 232547, 309527,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 386657, 0, 3,
                                                                       375317, 302075, 376937,
                                                                       232547, 233807, 311147,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 388682, 0, 3,
                                                                       378557, 304667, 380582,
                                                                       236327, 237867, 312767,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 391157, 0, 3,
                                                                       380582, 306287, 382607,
                                                                       237867, 239407, 314747,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 393632, 0, 3,
                                                                       382607, 307907, 384632,
                                                                       239407, 240947, 316727,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 396107, 0, 3,
                                                                       384632, 309527, 386657,
                                                                       240947, 242487, 318707,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 398582, 0, 3,
                                                                       388682, 312767, 391157,
                                                                       245567, 247415, 320687,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 401552, 0, 3,
                                                                       391157, 314747, 393632,
                                                                       247415, 249263, 323063,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 404522, 0, 3,
                                                                       393632, 316727, 396107,
                                                                       249263, 251111, 325439,
                                                                       ncols, gamma, p, q);

                    compute_prim_osl_three_center_electron_repulsion_0(buffer, 407492, 0, 3,
                                                                       398582, 320687, 401552,
                                                                       254807, 256991, 327815,
                                                                       ncols, gamma, p, q);

                    compute_prim_osl_three_center_electron_repulsion_0(buffer, 411002, 0, 3,
                                                                       401552, 323063, 404522,
                                                                       256991, 259175, 330623,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsl_three_center_electron_repulsion_0(buffer, 414512, 0, 3,
                                                                       407492, 327815, 411002,
                                                                       263543, 266091, 333431,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 418607, 360017, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 420343, 368837, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 422575, 378557, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 425365, 388682, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 428775, 398582, 2970, ncols);

                    simdfunc::contract_primitives(buffer, 432867, 407492, 3510, ncols);

                    simdfunc::contract_primitives(buffer, 437703, 414512, 4095, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 419867, 418607, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 421963, 420343, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 424600, 422575, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 427840, 425365, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 431745, 428775, 66, 1, nmax);

        simdtrf::transform_l_inner(buffer, 436377, 432867, 78, 1, nmax);

        simdtrf::transform_l_inner(buffer, 441798, 437703, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 443345, 419867, 421963, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 444773, 421963, 424600, 17,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 446609, 424600, 427840, 17,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 448904, 427840, 431745, 17,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 451709, 431745, 436377, 17,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 455075, 436377, 441798, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 459053, 443345, 444773, 17,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 461909, 444773, 446609, 17,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 465581, 446609, 448904, 17,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 470171, 448904, 451709, 17,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 475781, 451709, 455075, 17,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 482513, 459053, 461909, 17,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 487273, 461909, 465581, 17,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 493393, 465581, 470171, 17,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 501043, 470171, 475781, 17,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 510393, 482513, 487273, 17,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 517533, 487273, 493393, 17,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 526713, 493393, 501043, 17,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 538188, 510393, 517533, 17,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 548184, 517533, 526713, 17,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 561036, 538188, 548184, 17,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 574364, 561036, 28, 17, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 574364, 221, nmax);
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
