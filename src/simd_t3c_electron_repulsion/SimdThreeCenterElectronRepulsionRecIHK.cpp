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


#include "SimdThreeCenterElectronRepulsionRecIHK.hpp"

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
compute_ihk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ihk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 307055, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2145 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 307055, 222047, 14538, dimensions);

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
                                                        16, 17, 18}, ncols, fj, mu, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 7, 8,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 8, 9,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 9, 10,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 10, 11,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 11, 12,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 12, 13,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 13, 14,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 14, 15,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 124, 0, 3, 15, 16,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 130, 0, 3, 16, 17,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 136, 0, 3, 17, 18,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 18, 19,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 19, 20,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 154, 0, 3, 20, 21,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 160, 0, 3, 21, 22,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 166, 0, 3, 22, 23,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 25, 28,
                                                                       76, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 28, 31,
                                                                       82, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 31, 34,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 34, 37,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 37, 40,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 40, 43,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 43, 46,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 46, 49,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 252, 0, 3, 49, 52,
                                                                       124, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 262, 0, 3, 52, 55,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 55, 58,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 282, 0, 3, 58, 61,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 292, 0, 3, 61, 64,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 302, 0, 3, 64, 67,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 312, 0, 3, 67, 70,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 322, 0, 3, 76, 82,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 337, 0, 3, 82, 88,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 352, 0, 3, 88, 94,
                                                                       192, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 367, 0, 3, 94,
                                                                       100, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 382, 0, 3, 100,
                                                                       106, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 397, 0, 3, 106,
                                                                       112, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 412, 0, 3, 112,
                                                                       118, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 427, 0, 3, 118,
                                                                       124, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 442, 0, 3, 124,
                                                                       130, 252, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 457, 0, 3, 130,
                                                                       136, 262, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 472, 0, 3, 136,
                                                                       142, 272, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 487, 0, 3, 142,
                                                                       148, 282, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 502, 0, 3, 148,
                                                                       154, 292, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 517, 0, 3, 154,
                                                                       160, 302, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 532, 0, 3, 172,
                                                                       182, 322, 337, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 553, 0, 3, 182,
                                                                       192, 337, 352, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 574, 0, 3, 192,
                                                                       202, 352, 367, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 595, 0, 3, 202,
                                                                       212, 367, 382, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 616, 0, 3, 212,
                                                                       222, 382, 397, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 637, 0, 3, 222,
                                                                       232, 397, 412, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 658, 0, 3, 232,
                                                                       242, 412, 427, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 679, 0, 3, 242,
                                                                       252, 427, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 700, 0, 3, 252,
                                                                       262, 442, 457, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 721, 0, 3, 262,
                                                                       272, 457, 472, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 742, 0, 3, 272,
                                                                       282, 472, 487, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 763, 0, 3, 282,
                                                                       292, 487, 502, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 784, 0, 3, 292,
                                                                       302, 502, 517, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 805, 0, 3, 322,
                                                                       337, 532, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 833, 0, 3, 337,
                                                                       352, 553, 574, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 861, 0, 3, 352,
                                                                       367, 574, 595, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 889, 0, 3, 367,
                                                                       382, 595, 616, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 917, 0, 3, 382,
                                                                       397, 616, 637, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 945, 0, 3, 397,
                                                                       412, 637, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 973, 0, 3, 412,
                                                                       427, 658, 679, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1001, 0, 3, 427,
                                                                       442, 679, 700, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1029, 0, 3, 442,
                                                                       457, 700, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1057, 0, 3, 457,
                                                                       472, 721, 742, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1085, 0, 3, 472,
                                                                       487, 742, 763, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1113, 0, 3, 487,
                                                                       502, 763, 784, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1141, 0, 3, 532,
                                                                       553, 805, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1177, 0, 3, 553,
                                                                       574, 833, 861, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1213, 0, 3, 574,
                                                                       595, 861, 889, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1249, 0, 3, 595,
                                                                       616, 889, 917, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1285, 0, 3, 616,
                                                                       637, 917, 945, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1321, 0, 3, 637,
                                                                       658, 945, 973, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1357, 0, 3, 658,
                                                                       679, 973, 1001, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1393, 0, 3, 679,
                                                                       700, 1001, 1029, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1429, 0, 3, 700,
                                                                       721, 1029, 1057, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1465, 0, 3, 721,
                                                                       742, 1057, 1085, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1501, 0, 3, 742,
                                                                       763, 1085, 1113, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1537, 0, 3, 805,
                                                                       833, 1141, 1177, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1582, 0, 3, 833,
                                                                       861, 1177, 1213, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1627, 0, 3, 861,
                                                                       889, 1213, 1249, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1672, 0, 3, 889,
                                                                       917, 1249, 1285, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1717, 0, 3, 917,
                                                                       945, 1285, 1321, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1762, 0, 3, 945,
                                                                       973, 1321, 1357, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1807, 0, 3, 973,
                                                                       1001, 1357, 1393, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1852, 0, 3, 1001,
                                                                       1029, 1393, 1429, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1897, 0, 3, 1029,
                                                                       1057, 1429, 1465, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1942, 0, 3, 1057,
                                                                       1085, 1465, 1501, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1987, 0, 3, 1141,
                                                                       1177, 1537, 1582, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2042, 0, 3, 1177,
                                                                       1213, 1582, 1627, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2097, 0, 3, 1213,
                                                                       1249, 1627, 1672, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2152, 0, 3, 1249,
                                                                       1285, 1672, 1717, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2207, 0, 3, 1285,
                                                                       1321, 1717, 1762, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2262, 0, 3, 1321,
                                                                       1357, 1762, 1807, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2317, 0, 3, 1357,
                                                                       1393, 1807, 1852, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1393,
                                                                       1429, 1852, 1897, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2427, 0, 3, 1429,
                                                                       1465, 1897, 1942, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2482, 0, 3, 1537,
                                                                       1582, 1987, 2042, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2548, 0, 3, 1582,
                                                                       1627, 2042, 2097, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2614, 0, 3, 1627,
                                                                       1672, 2097, 2152, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2680, 0, 3, 1672,
                                                                       1717, 2152, 2207, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2746, 0, 3, 1717,
                                                                       1762, 2207, 2262, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2812, 0, 3, 1762,
                                                                       1807, 2262, 2317, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2878, 0, 3, 1807,
                                                                       1852, 2317, 2372, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2944, 0, 3, 1852,
                                                                       1897, 2372, 2427, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3010, 0, 3, 1987,
                                                                       2042, 2482, 2548, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 2042,
                                                                       2097, 2548, 2614, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3166, 0, 3, 2097,
                                                                       2152, 2614, 2680, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3244, 0, 3, 2152,
                                                                       2207, 2680, 2746, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3322, 0, 3, 2207,
                                                                       2262, 2746, 2812, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3400, 0, 3, 2262,
                                                                       2317, 2812, 2878, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3478, 0, 3, 2317,
                                                                       2372, 2878, 2944, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3556, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3559, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3562, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3565, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3568, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3571, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3574, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3577, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3580, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3583, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3586, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3589, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3592, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3595, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3598, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3601, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3604, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3607, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3610, 3, 7, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3619, 3, 8, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3628, 3, 9, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3637, 3, 10, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3646, 3, 11, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3655, 3, 12, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3664, 3, 13, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3673, 3, 14, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3682, 3, 15, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3691, 3, 16, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3700, 3, 17, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3709, 3, 18, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3718, 3, 19, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3727, 3, 20, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3736, 3, 21, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3745, 3, 22, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3754, 3, 23, 73,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3763, 3, 25, 76,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3781, 3, 28, 82,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3799, 3, 31, 88,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3817, 3, 34, 94,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3835, 3, 37, 100,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3853, 3, 40, 106,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3871, 3, 43, 112,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3889, 3, 46, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3907, 3, 49, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3925, 3, 52, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3943, 3, 55, 136,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3961, 3, 58, 142,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3979, 3, 61, 148,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3997, 3, 64, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4015, 3, 67, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4033, 3, 70, 166,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4051, 3, 76, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4081, 3, 82, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4111, 3, 88, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4141, 3, 94, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4171, 3, 100, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4201, 3, 106, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4231, 3, 112, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4261, 3, 118, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4291, 3, 124, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4321, 3, 130, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4351, 3, 136, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4381, 3, 142, 282,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4411, 3, 148, 292,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4441, 3, 154, 302,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4471, 3, 160, 312,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4501, 3, 172, 322,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4546, 3, 182, 337,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4591, 3, 192, 352,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4636, 3, 202, 367,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4681, 3, 212, 382,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4726, 3, 222, 397,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4771, 3, 232, 412,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4816, 3, 242, 427,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4861, 3, 252, 442,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4906, 3, 262, 457,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4951, 3, 272, 472,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4996, 3, 282, 487,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5041, 3, 292, 502,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5086, 3, 302, 517,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5131, 3, 322, 532,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5194, 3, 337, 553,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5257, 3, 352, 574,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5320, 3, 367, 595,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5383, 3, 382, 616,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5446, 3, 397, 637,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5509, 3, 412, 658,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5572, 3, 427, 679,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5635, 3, 442, 700,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5698, 3, 457, 721,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5761, 3, 472, 742,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5824, 3, 487, 763,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5887, 3, 502, 784,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5950, 3, 532, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6034, 3, 553, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6118, 3, 574, 861,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6202, 3, 595, 889,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6286, 3, 616, 917,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6370, 3, 637, 945,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6454, 3, 658, 973,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6538, 3, 679,
                                                                       1001, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6622, 3, 700,
                                                                       1029, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6706, 3, 721,
                                                                       1057, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6790, 3, 742,
                                                                       1085, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6874, 3, 763,
                                                                       1113, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6958, 3, 805,
                                                                       1141, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7066, 3, 833,
                                                                       1177, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7174, 3, 861,
                                                                       1213, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7282, 3, 889,
                                                                       1249, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7390, 3, 917,
                                                                       1285, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7498, 3, 945,
                                                                       1321, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7606, 3, 973,
                                                                       1357, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7714, 3, 1001,
                                                                       1393, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7822, 3, 1029,
                                                                       1429, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7930, 3, 1057,
                                                                       1465, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8038, 3, 1085,
                                                                       1501, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8146, 3, 1141,
                                                                       1537, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8281, 3, 1177,
                                                                       1582, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8416, 3, 1213,
                                                                       1627, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8551, 3, 1249,
                                                                       1672, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8686, 3, 1285,
                                                                       1717, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8821, 3, 1321,
                                                                       1762, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8956, 3, 1357,
                                                                       1807, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9091, 3, 1393,
                                                                       1852, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9226, 3, 1429,
                                                                       1897, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9361, 3, 1465,
                                                                       1942, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9496, 3, 1537,
                                                                       1987, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9661, 3, 1582,
                                                                       2042, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9826, 3, 1627,
                                                                       2097, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9991, 3, 1672,
                                                                       2152, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10156, 3, 1717,
                                                                       2207, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10321, 3, 1762,
                                                                       2262, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10486, 3, 1807,
                                                                       2317, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10651, 3, 1852,
                                                                       2372, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10816, 3, 1897,
                                                                       2427, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10981, 3, 1987,
                                                                       2482, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11179, 3, 2042,
                                                                       2548, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11377, 3, 2097,
                                                                       2614, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11575, 3, 2152,
                                                                       2680, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11773, 3, 2207,
                                                                       2746, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11971, 3, 2262,
                                                                       2812, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12169, 3, 2317,
                                                                       2878, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12367, 3, 2372,
                                                                       2944, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12565, 3, 2482,
                                                                       3010, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12799, 3, 2548,
                                                                       3088, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13033, 3, 2614,
                                                                       3166, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13267, 3, 2680,
                                                                       3244, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13501, 3, 2746,
                                                                       3322, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13735, 3, 2812,
                                                                       3400, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13969, 3, 2878,
                                                                       3478, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14203, 3, 7, 8,
                                                                       3562, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14209, 3, 8, 9,
                                                                       3565, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14215, 3, 9, 10,
                                                                       3568, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14221, 3, 10, 11,
                                                                       3571, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14227, 3, 11, 12,
                                                                       3574, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14233, 3, 12, 13,
                                                                       3577, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14239, 3, 13, 14,
                                                                       3580, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14245, 3, 14, 15,
                                                                       3583, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14251, 3, 15, 16,
                                                                       3586, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14257, 3, 16, 17,
                                                                       3589, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14263, 3, 17, 18,
                                                                       3592, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14269, 3, 18, 19,
                                                                       3595, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14275, 3, 19, 20,
                                                                       3598, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14281, 3, 20, 21,
                                                                       3601, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14287, 3, 21, 22,
                                                                       3604, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14293, 3, 22, 23,
                                                                       3607, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14299, 0, 3,
                                                                       14203, 3562, 14209, 3628,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14317, 0, 3,
                                                                       14209, 3565, 14215, 3637,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14335, 0, 3,
                                                                       14215, 3568, 14221, 3646,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14353, 0, 3,
                                                                       14221, 3571, 14227, 3655,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14371, 0, 3,
                                                                       14227, 3574, 14233, 3664,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14389, 0, 3,
                                                                       14233, 3577, 14239, 3673,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14407, 0, 3,
                                                                       14239, 3580, 14245, 3682,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14425, 0, 3,
                                                                       14245, 3583, 14251, 3691,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14443, 0, 3,
                                                                       14251, 3586, 14257, 3700,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14461, 0, 3,
                                                                       14257, 3589, 14263, 3709,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14479, 0, 3,
                                                                       14263, 3592, 14269, 3718,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14497, 0, 3,
                                                                       14269, 3595, 14275, 3727,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14515, 0, 3,
                                                                       14275, 3598, 14281, 3736,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14533, 0, 3,
                                                                       14281, 3601, 14287, 3745,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14551, 0, 3,
                                                                       14287, 3604, 14293, 3754,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14569, 0, 3,
                                                                       14299, 3628, 14317, 76,
                                                                       82, 3799, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14605, 0, 3,
                                                                       14317, 3637, 14335, 82,
                                                                       88, 3817, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14641, 0, 3,
                                                                       14335, 3646, 14353, 88,
                                                                       94, 3835, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14677, 0, 3,
                                                                       14353, 3655, 14371, 94,
                                                                       100, 3853, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14713, 0, 3,
                                                                       14371, 3664, 14389, 100,
                                                                       106, 3871, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14749, 0, 3,
                                                                       14389, 3673, 14407, 106,
                                                                       112, 3889, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14785, 0, 3,
                                                                       14407, 3682, 14425, 112,
                                                                       118, 3907, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14821, 0, 3,
                                                                       14425, 3691, 14443, 118,
                                                                       124, 3925, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14857, 0, 3,
                                                                       14443, 3700, 14461, 124,
                                                                       130, 3943, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14893, 0, 3,
                                                                       14461, 3709, 14479, 130,
                                                                       136, 3961, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14929, 0, 3,
                                                                       14479, 3718, 14497, 136,
                                                                       142, 3979, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14965, 0, 3,
                                                                       14497, 3727, 14515, 142,
                                                                       148, 3997, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15001, 0, 3,
                                                                       14515, 3736, 14533, 148,
                                                                       154, 4015, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15037, 0, 3,
                                                                       14533, 3745, 14551, 154,
                                                                       160, 4033, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15073, 0, 3,
                                                                       14569, 3799, 14605, 172,
                                                                       182, 4111, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15133, 0, 3,
                                                                       14605, 3817, 14641, 182,
                                                                       192, 4141, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15193, 0, 3,
                                                                       14641, 3835, 14677, 192,
                                                                       202, 4171, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15253, 0, 3,
                                                                       14677, 3853, 14713, 202,
                                                                       212, 4201, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15313, 0, 3,
                                                                       14713, 3871, 14749, 212,
                                                                       222, 4231, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15373, 0, 3,
                                                                       14749, 3889, 14785, 222,
                                                                       232, 4261, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15433, 0, 3,
                                                                       14785, 3907, 14821, 232,
                                                                       242, 4291, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15493, 0, 3,
                                                                       14821, 3925, 14857, 242,
                                                                       252, 4321, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15553, 0, 3,
                                                                       14857, 3943, 14893, 252,
                                                                       262, 4351, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15613, 0, 3,
                                                                       14893, 3961, 14929, 262,
                                                                       272, 4381, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15673, 0, 3,
                                                                       14929, 3979, 14965, 272,
                                                                       282, 4411, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15733, 0, 3,
                                                                       14965, 3997, 15001, 282,
                                                                       292, 4441, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15793, 0, 3,
                                                                       15001, 4015, 15037, 292,
                                                                       302, 4471, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15853, 0, 3,
                                                                       15073, 4111, 15133, 322,
                                                                       337, 4591, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15943, 0, 3,
                                                                       15133, 4141, 15193, 337,
                                                                       352, 4636, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16033, 0, 3,
                                                                       15193, 4171, 15253, 352,
                                                                       367, 4681, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16123, 0, 3,
                                                                       15253, 4201, 15313, 367,
                                                                       382, 4726, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16213, 0, 3,
                                                                       15313, 4231, 15373, 382,
                                                                       397, 4771, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16303, 0, 3,
                                                                       15373, 4261, 15433, 397,
                                                                       412, 4816, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16393, 0, 3,
                                                                       15433, 4291, 15493, 412,
                                                                       427, 4861, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16483, 0, 3,
                                                                       15493, 4321, 15553, 427,
                                                                       442, 4906, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16573, 0, 3,
                                                                       15553, 4351, 15613, 442,
                                                                       457, 4951, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16663, 0, 3,
                                                                       15613, 4381, 15673, 457,
                                                                       472, 4996, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16753, 0, 3,
                                                                       15673, 4411, 15733, 472,
                                                                       487, 5041, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16843, 0, 3,
                                                                       15733, 4441, 15793, 487,
                                                                       502, 5086, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16933, 0, 3,
                                                                       15853, 4591, 15943, 532,
                                                                       553, 5257, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17059, 0, 3,
                                                                       15943, 4636, 16033, 553,
                                                                       574, 5320, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17185, 0, 3,
                                                                       16033, 4681, 16123, 574,
                                                                       595, 5383, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17311, 0, 3,
                                                                       16123, 4726, 16213, 595,
                                                                       616, 5446, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17437, 0, 3,
                                                                       16213, 4771, 16303, 616,
                                                                       637, 5509, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17563, 0, 3,
                                                                       16303, 4816, 16393, 637,
                                                                       658, 5572, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17689, 0, 3,
                                                                       16393, 4861, 16483, 658,
                                                                       679, 5635, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17815, 0, 3,
                                                                       16483, 4906, 16573, 679,
                                                                       700, 5698, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17941, 0, 3,
                                                                       16573, 4951, 16663, 700,
                                                                       721, 5761, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18067, 0, 3,
                                                                       16663, 4996, 16753, 721,
                                                                       742, 5824, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18193, 0, 3,
                                                                       16753, 5041, 16843, 742,
                                                                       763, 5887, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18319, 0, 3,
                                                                       16933, 5257, 17059, 805,
                                                                       833, 6118, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18487, 0, 3,
                                                                       17059, 5320, 17185, 833,
                                                                       861, 6202, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18655, 0, 3,
                                                                       17185, 5383, 17311, 861,
                                                                       889, 6286, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18823, 0, 3,
                                                                       17311, 5446, 17437, 889,
                                                                       917, 6370, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18991, 0, 3,
                                                                       17437, 5509, 17563, 917,
                                                                       945, 6454, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19159, 0, 3,
                                                                       17563, 5572, 17689, 945,
                                                                       973, 6538, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19327, 0, 3,
                                                                       17689, 5635, 17815, 973,
                                                                       1001, 6622, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19495, 0, 3,
                                                                       17815, 5698, 17941, 1001,
                                                                       1029, 6706, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19663, 0, 3,
                                                                       17941, 5761, 18067, 1029,
                                                                       1057, 6790, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19831, 0, 3,
                                                                       18067, 5824, 18193, 1057,
                                                                       1085, 6874, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19999, 0, 3,
                                                                       18319, 6118, 18487, 1141,
                                                                       1177, 7174, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20215, 0, 3,
                                                                       18487, 6202, 18655, 1177,
                                                                       1213, 7282, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20431, 0, 3,
                                                                       18655, 6286, 18823, 1213,
                                                                       1249, 7390, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20647, 0, 3,
                                                                       18823, 6370, 18991, 1249,
                                                                       1285, 7498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20863, 0, 3,
                                                                       18991, 6454, 19159, 1285,
                                                                       1321, 7606, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21079, 0, 3,
                                                                       19159, 6538, 19327, 1321,
                                                                       1357, 7714, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21295, 0, 3,
                                                                       19327, 6622, 19495, 1357,
                                                                       1393, 7822, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21511, 0, 3,
                                                                       19495, 6706, 19663, 1393,
                                                                       1429, 7930, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21727, 0, 3,
                                                                       19663, 6790, 19831, 1429,
                                                                       1465, 8038, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21943, 0, 3,
                                                                       19999, 7174, 20215, 1537,
                                                                       1582, 8416, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22213, 0, 3,
                                                                       20215, 7282, 20431, 1582,
                                                                       1627, 8551, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22483, 0, 3,
                                                                       20431, 7390, 20647, 1627,
                                                                       1672, 8686, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22753, 0, 3,
                                                                       20647, 7498, 20863, 1672,
                                                                       1717, 8821, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23023, 0, 3,
                                                                       20863, 7606, 21079, 1717,
                                                                       1762, 8956, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23293, 0, 3,
                                                                       21079, 7714, 21295, 1762,
                                                                       1807, 9091, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23563, 0, 3,
                                                                       21295, 7822, 21511, 1807,
                                                                       1852, 9226, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23833, 0, 3,
                                                                       21511, 7930, 21727, 1852,
                                                                       1897, 9361, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24103, 0, 3,
                                                                       21943, 8416, 22213, 1987,
                                                                       2042, 9826, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24433, 0, 3,
                                                                       22213, 8551, 22483, 2042,
                                                                       2097, 9991, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24763, 0, 3,
                                                                       22483, 8686, 22753, 2097,
                                                                       2152, 10156, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25093, 0, 3,
                                                                       22753, 8821, 23023, 2152,
                                                                       2207, 10321, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25423, 0, 3,
                                                                       23023, 8956, 23293, 2207,
                                                                       2262, 10486, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25753, 0, 3,
                                                                       23293, 9091, 23563, 2262,
                                                                       2317, 10651, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 26083, 0, 3,
                                                                       23563, 9226, 23833, 2317,
                                                                       2372, 10816, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26413, 0, 3,
                                                                       24103, 9826, 24433, 2482,
                                                                       2548, 11377, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26809, 0, 3,
                                                                       24433, 9991, 24763, 2548,
                                                                       2614, 11575, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27205, 0, 3,
                                                                       24763, 10156, 25093, 2614,
                                                                       2680, 11773, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27601, 0, 3,
                                                                       25093, 10321, 25423, 2680,
                                                                       2746, 11971, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27997, 0, 3,
                                                                       25423, 10486, 25753, 2746,
                                                                       2812, 12169, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 28393, 0, 3,
                                                                       25753, 10651, 26083, 2812,
                                                                       2878, 12367, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 28789, 0, 3,
                                                                       26413, 11377, 26809, 3010,
                                                                       3088, 13033, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 29257, 0, 3,
                                                                       26809, 11575, 27205, 3088,
                                                                       3166, 13267, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 29725, 0, 3,
                                                                       27205, 11773, 27601, 3166,
                                                                       3244, 13501, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 30193, 0, 3,
                                                                       27601, 11971, 27997, 3244,
                                                                       3322, 13735, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 30661, 0, 3,
                                                                       27997, 12169, 28393, 3322,
                                                                       3400, 13969, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31129, 3, 3556,
                                                                       3559, 14203, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31139, 3, 3559,
                                                                       3562, 14209, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31149, 3, 3562,
                                                                       3565, 14215, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31159, 3, 3565,
                                                                       3568, 14221, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31169, 3, 3568,
                                                                       3571, 14227, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31179, 3, 3571,
                                                                       3574, 14233, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31189, 3, 3574,
                                                                       3577, 14239, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31199, 3, 3577,
                                                                       3580, 14245, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31209, 3, 3580,
                                                                       3583, 14251, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31219, 3, 3583,
                                                                       3586, 14257, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31229, 3, 3586,
                                                                       3589, 14263, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31239, 3, 3589,
                                                                       3592, 14269, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31249, 3, 3592,
                                                                       3595, 14275, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31259, 3, 3595,
                                                                       3598, 14281, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31269, 3, 3598,
                                                                       3601, 14287, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31279, 3, 3601,
                                                                       3604, 14293, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31289, 0, 3,
                                                                       31129, 14203, 31139, 3610,
                                                                       3619, 14299, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31319, 0, 3,
                                                                       31139, 14209, 31149, 3619,
                                                                       3628, 14317, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31349, 0, 3,
                                                                       31149, 14215, 31159, 3628,
                                                                       3637, 14335, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31379, 0, 3,
                                                                       31159, 14221, 31169, 3637,
                                                                       3646, 14353, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31409, 0, 3,
                                                                       31169, 14227, 31179, 3646,
                                                                       3655, 14371, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31439, 0, 3,
                                                                       31179, 14233, 31189, 3655,
                                                                       3664, 14389, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31469, 0, 3,
                                                                       31189, 14239, 31199, 3664,
                                                                       3673, 14407, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31499, 0, 3,
                                                                       31199, 14245, 31209, 3673,
                                                                       3682, 14425, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31529, 0, 3,
                                                                       31209, 14251, 31219, 3682,
                                                                       3691, 14443, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31559, 0, 3,
                                                                       31219, 14257, 31229, 3691,
                                                                       3700, 14461, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31589, 0, 3,
                                                                       31229, 14263, 31239, 3700,
                                                                       3709, 14479, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31619, 0, 3,
                                                                       31239, 14269, 31249, 3709,
                                                                       3718, 14497, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31649, 0, 3,
                                                                       31249, 14275, 31259, 3718,
                                                                       3727, 14515, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31679, 0, 3,
                                                                       31259, 14281, 31269, 3727,
                                                                       3736, 14533, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31709, 0, 3,
                                                                       31269, 14287, 31279, 3736,
                                                                       3745, 14551, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31739, 0, 3,
                                                                       31289, 14299, 31319, 3763,
                                                                       3781, 14569, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31799, 0, 3,
                                                                       31319, 14317, 31349, 3781,
                                                                       3799, 14605, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31859, 0, 3,
                                                                       31349, 14335, 31379, 3799,
                                                                       3817, 14641, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31919, 0, 3,
                                                                       31379, 14353, 31409, 3817,
                                                                       3835, 14677, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31979, 0, 3,
                                                                       31409, 14371, 31439, 3835,
                                                                       3853, 14713, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32039, 0, 3,
                                                                       31439, 14389, 31469, 3853,
                                                                       3871, 14749, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32099, 0, 3,
                                                                       31469, 14407, 31499, 3871,
                                                                       3889, 14785, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32159, 0, 3,
                                                                       31499, 14425, 31529, 3889,
                                                                       3907, 14821, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32219, 0, 3,
                                                                       31529, 14443, 31559, 3907,
                                                                       3925, 14857, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32279, 0, 3,
                                                                       31559, 14461, 31589, 3925,
                                                                       3943, 14893, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32339, 0, 3,
                                                                       31589, 14479, 31619, 3943,
                                                                       3961, 14929, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32399, 0, 3,
                                                                       31619, 14497, 31649, 3961,
                                                                       3979, 14965, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32459, 0, 3,
                                                                       31649, 14515, 31679, 3979,
                                                                       3997, 15001, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32519, 0, 3,
                                                                       31679, 14533, 31709, 3997,
                                                                       4015, 15037, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32579, 0, 3,
                                                                       31739, 14569, 31799, 4051,
                                                                       4081, 15073, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32679, 0, 3,
                                                                       31799, 14605, 31859, 4081,
                                                                       4111, 15133, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32779, 0, 3,
                                                                       31859, 14641, 31919, 4111,
                                                                       4141, 15193, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32879, 0, 3,
                                                                       31919, 14677, 31979, 4141,
                                                                       4171, 15253, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32979, 0, 3,
                                                                       31979, 14713, 32039, 4171,
                                                                       4201, 15313, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33079, 0, 3,
                                                                       32039, 14749, 32099, 4201,
                                                                       4231, 15373, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33179, 0, 3,
                                                                       32099, 14785, 32159, 4231,
                                                                       4261, 15433, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33279, 0, 3,
                                                                       32159, 14821, 32219, 4261,
                                                                       4291, 15493, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33379, 0, 3,
                                                                       32219, 14857, 32279, 4291,
                                                                       4321, 15553, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33479, 0, 3,
                                                                       32279, 14893, 32339, 4321,
                                                                       4351, 15613, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33579, 0, 3,
                                                                       32339, 14929, 32399, 4351,
                                                                       4381, 15673, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33679, 0, 3,
                                                                       32399, 14965, 32459, 4381,
                                                                       4411, 15733, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33779, 0, 3,
                                                                       32459, 15001, 32519, 4411,
                                                                       4441, 15793, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 33879, 0, 3,
                                                                       32579, 15073, 32679, 4501,
                                                                       4546, 15853, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34029, 0, 3,
                                                                       32679, 15133, 32779, 4546,
                                                                       4591, 15943, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34179, 0, 3,
                                                                       32779, 15193, 32879, 4591,
                                                                       4636, 16033, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34329, 0, 3,
                                                                       32879, 15253, 32979, 4636,
                                                                       4681, 16123, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34479, 0, 3,
                                                                       32979, 15313, 33079, 4681,
                                                                       4726, 16213, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34629, 0, 3,
                                                                       33079, 15373, 33179, 4726,
                                                                       4771, 16303, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34779, 0, 3,
                                                                       33179, 15433, 33279, 4771,
                                                                       4816, 16393, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34929, 0, 3,
                                                                       33279, 15493, 33379, 4816,
                                                                       4861, 16483, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35079, 0, 3,
                                                                       33379, 15553, 33479, 4861,
                                                                       4906, 16573, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35229, 0, 3,
                                                                       33479, 15613, 33579, 4906,
                                                                       4951, 16663, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35379, 0, 3,
                                                                       33579, 15673, 33679, 4951,
                                                                       4996, 16753, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35529, 0, 3,
                                                                       33679, 15733, 33779, 4996,
                                                                       5041, 16843, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35679, 0, 3,
                                                                       33879, 15853, 34029, 5131,
                                                                       5194, 16933, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35889, 0, 3,
                                                                       34029, 15943, 34179, 5194,
                                                                       5257, 17059, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36099, 0, 3,
                                                                       34179, 16033, 34329, 5257,
                                                                       5320, 17185, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36309, 0, 3,
                                                                       34329, 16123, 34479, 5320,
                                                                       5383, 17311, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36519, 0, 3,
                                                                       34479, 16213, 34629, 5383,
                                                                       5446, 17437, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36729, 0, 3,
                                                                       34629, 16303, 34779, 5446,
                                                                       5509, 17563, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36939, 0, 3,
                                                                       34779, 16393, 34929, 5509,
                                                                       5572, 17689, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37149, 0, 3,
                                                                       34929, 16483, 35079, 5572,
                                                                       5635, 17815, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37359, 0, 3,
                                                                       35079, 16573, 35229, 5635,
                                                                       5698, 17941, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37569, 0, 3,
                                                                       35229, 16663, 35379, 5698,
                                                                       5761, 18067, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37779, 0, 3,
                                                                       35379, 16753, 35529, 5761,
                                                                       5824, 18193, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37989, 0, 3,
                                                                       35679, 16933, 35889, 5950,
                                                                       6034, 18319, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38269, 0, 3,
                                                                       35889, 17059, 36099, 6034,
                                                                       6118, 18487, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38549, 0, 3,
                                                                       36099, 17185, 36309, 6118,
                                                                       6202, 18655, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38829, 0, 3,
                                                                       36309, 17311, 36519, 6202,
                                                                       6286, 18823, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39109, 0, 3,
                                                                       36519, 17437, 36729, 6286,
                                                                       6370, 18991, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39389, 0, 3,
                                                                       36729, 17563, 36939, 6370,
                                                                       6454, 19159, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39669, 0, 3,
                                                                       36939, 17689, 37149, 6454,
                                                                       6538, 19327, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39949, 0, 3,
                                                                       37149, 17815, 37359, 6538,
                                                                       6622, 19495, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 40229, 0, 3,
                                                                       37359, 17941, 37569, 6622,
                                                                       6706, 19663, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 40509, 0, 3,
                                                                       37569, 18067, 37779, 6706,
                                                                       6790, 19831, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40789, 0, 3,
                                                                       37989, 18319, 38269, 6958,
                                                                       7066, 19999, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41149, 0, 3,
                                                                       38269, 18487, 38549, 7066,
                                                                       7174, 20215, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41509, 0, 3,
                                                                       38549, 18655, 38829, 7174,
                                                                       7282, 20431, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41869, 0, 3,
                                                                       38829, 18823, 39109, 7282,
                                                                       7390, 20647, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42229, 0, 3,
                                                                       39109, 18991, 39389, 7390,
                                                                       7498, 20863, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42589, 0, 3,
                                                                       39389, 19159, 39669, 7498,
                                                                       7606, 21079, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42949, 0, 3,
                                                                       39669, 19327, 39949, 7606,
                                                                       7714, 21295, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 43309, 0, 3,
                                                                       39949, 19495, 40229, 7714,
                                                                       7822, 21511, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 43669, 0, 3,
                                                                       40229, 19663, 40509, 7822,
                                                                       7930, 21727, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44029, 0, 3,
                                                                       40789, 19999, 41149, 8146,
                                                                       8281, 21943, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44479, 0, 3,
                                                                       41149, 20215, 41509, 8281,
                                                                       8416, 22213, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44929, 0, 3,
                                                                       41509, 20431, 41869, 8416,
                                                                       8551, 22483, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 45379, 0, 3,
                                                                       41869, 20647, 42229, 8551,
                                                                       8686, 22753, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 45829, 0, 3,
                                                                       42229, 20863, 42589, 8686,
                                                                       8821, 23023, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 46279, 0, 3,
                                                                       42589, 21079, 42949, 8821,
                                                                       8956, 23293, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 46729, 0, 3,
                                                                       42949, 21295, 43309, 8956,
                                                                       9091, 23563, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 47179, 0, 3,
                                                                       43309, 21511, 43669, 9091,
                                                                       9226, 23833, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 47629, 0, 3,
                                                                       44029, 21943, 44479, 9496,
                                                                       9661, 24103, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 48179, 0, 3,
                                                                       44479, 22213, 44929, 9661,
                                                                       9826, 24433, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 48729, 0, 3,
                                                                       44929, 22483, 45379, 9826,
                                                                       9991, 24763, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 49279, 0, 3,
                                                                       45379, 22753, 45829, 9991,
                                                                       10156, 25093, ncols,
                                                                       gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 49829, 0, 3,
                                                                       45829, 23023, 46279,
                                                                       10156, 10321, 25423,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 50379, 0, 3,
                                                                       46279, 23293, 46729,
                                                                       10321, 10486, 25753,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 50929, 0, 3,
                                                                       46729, 23563, 47179,
                                                                       10486, 10651, 26083,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 51479, 0, 3,
                                                                       47629, 24103, 48179,
                                                                       10981, 11179, 26413,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 52139, 0, 3,
                                                                       48179, 24433, 48729,
                                                                       11179, 11377, 26809,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 52799, 0, 3,
                                                                       48729, 24763, 49279,
                                                                       11377, 11575, 27205,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 53459, 0, 3,
                                                                       49279, 25093, 49829,
                                                                       11575, 11773, 27601,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 54119, 0, 3,
                                                                       49829, 25423, 50379,
                                                                       11773, 11971, 27997,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 54779, 0, 3,
                                                                       50379, 25753, 50929,
                                                                       11971, 12169, 28393,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 55439, 0, 3,
                                                                       51479, 26413, 52139,
                                                                       12565, 12799, 28789,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 56219, 0, 3,
                                                                       52139, 26809, 52799,
                                                                       12799, 13033, 29257,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 56999, 0, 3,
                                                                       52799, 27205, 53459,
                                                                       13033, 13267, 29725,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 57779, 0, 3,
                                                                       53459, 27601, 54119,
                                                                       13267, 13501, 30193,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 58559, 0, 3,
                                                                       54119, 27997, 54779,
                                                                       13501, 13735, 30661,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59339, 3, 14203,
                                                                       14209, 31149, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59354, 3, 14209,
                                                                       14215, 31159, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59369, 3, 14215,
                                                                       14221, 31169, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59384, 3, 14221,
                                                                       14227, 31179, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59399, 3, 14227,
                                                                       14233, 31189, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59414, 3, 14233,
                                                                       14239, 31199, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59429, 3, 14239,
                                                                       14245, 31209, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59444, 3, 14245,
                                                                       14251, 31219, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59459, 3, 14251,
                                                                       14257, 31229, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59474, 3, 14257,
                                                                       14263, 31239, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59489, 3, 14263,
                                                                       14269, 31249, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59504, 3, 14269,
                                                                       14275, 31259, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59519, 3, 14275,
                                                                       14281, 31269, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59534, 3, 14281,
                                                                       14287, 31279, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59549, 0, 3,
                                                                       59339, 31149, 59354,
                                                                       14299, 14317, 31349,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59594, 0, 3,
                                                                       59354, 31159, 59369,
                                                                       14317, 14335, 31379,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59639, 0, 3,
                                                                       59369, 31169, 59384,
                                                                       14335, 14353, 31409,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59684, 0, 3,
                                                                       59384, 31179, 59399,
                                                                       14353, 14371, 31439,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59729, 0, 3,
                                                                       59399, 31189, 59414,
                                                                       14371, 14389, 31469,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59774, 0, 3,
                                                                       59414, 31199, 59429,
                                                                       14389, 14407, 31499,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59819, 0, 3,
                                                                       59429, 31209, 59444,
                                                                       14407, 14425, 31529,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59864, 0, 3,
                                                                       59444, 31219, 59459,
                                                                       14425, 14443, 31559,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59909, 0, 3,
                                                                       59459, 31229, 59474,
                                                                       14443, 14461, 31589,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59954, 0, 3,
                                                                       59474, 31239, 59489,
                                                                       14461, 14479, 31619,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59999, 0, 3,
                                                                       59489, 31249, 59504,
                                                                       14479, 14497, 31649,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 60044, 0, 3,
                                                                       59504, 31259, 59519,
                                                                       14497, 14515, 31679,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 60089, 0, 3,
                                                                       59519, 31269, 59534,
                                                                       14515, 14533, 31709,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60134, 0, 3,
                                                                       59549, 31349, 59594,
                                                                       14569, 14605, 31859,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60224, 0, 3,
                                                                       59594, 31379, 59639,
                                                                       14605, 14641, 31919,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60314, 0, 3,
                                                                       59639, 31409, 59684,
                                                                       14641, 14677, 31979,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60404, 0, 3,
                                                                       59684, 31439, 59729,
                                                                       14677, 14713, 32039,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60494, 0, 3,
                                                                       59729, 31469, 59774,
                                                                       14713, 14749, 32099,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60584, 0, 3,
                                                                       59774, 31499, 59819,
                                                                       14749, 14785, 32159,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60674, 0, 3,
                                                                       59819, 31529, 59864,
                                                                       14785, 14821, 32219,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60764, 0, 3,
                                                                       59864, 31559, 59909,
                                                                       14821, 14857, 32279,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60854, 0, 3,
                                                                       59909, 31589, 59954,
                                                                       14857, 14893, 32339,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60944, 0, 3,
                                                                       59954, 31619, 59999,
                                                                       14893, 14929, 32399,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 61034, 0, 3,
                                                                       59999, 31649, 60044,
                                                                       14929, 14965, 32459,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 61124, 0, 3,
                                                                       60044, 31679, 60089,
                                                                       14965, 15001, 32519,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61214, 0, 3,
                                                                       60134, 31859, 60224,
                                                                       15073, 15133, 32779,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61364, 0, 3,
                                                                       60224, 31919, 60314,
                                                                       15133, 15193, 32879,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61514, 0, 3,
                                                                       60314, 31979, 60404,
                                                                       15193, 15253, 32979,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61664, 0, 3,
                                                                       60404, 32039, 60494,
                                                                       15253, 15313, 33079,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61814, 0, 3,
                                                                       60494, 32099, 60584,
                                                                       15313, 15373, 33179,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61964, 0, 3,
                                                                       60584, 32159, 60674,
                                                                       15373, 15433, 33279,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 62114, 0, 3,
                                                                       60674, 32219, 60764,
                                                                       15433, 15493, 33379,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 62264, 0, 3,
                                                                       60764, 32279, 60854,
                                                                       15493, 15553, 33479,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 62414, 0, 3,
                                                                       60854, 32339, 60944,
                                                                       15553, 15613, 33579,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 62564, 0, 3,
                                                                       60944, 32399, 61034,
                                                                       15613, 15673, 33679,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 62714, 0, 3,
                                                                       61034, 32459, 61124,
                                                                       15673, 15733, 33779,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62864, 0, 3,
                                                                       61214, 32779, 61364,
                                                                       15853, 15943, 34179,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63089, 0, 3,
                                                                       61364, 32879, 61514,
                                                                       15943, 16033, 34329,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63314, 0, 3,
                                                                       61514, 32979, 61664,
                                                                       16033, 16123, 34479,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63539, 0, 3,
                                                                       61664, 33079, 61814,
                                                                       16123, 16213, 34629,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63764, 0, 3,
                                                                       61814, 33179, 61964,
                                                                       16213, 16303, 34779,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63989, 0, 3,
                                                                       61964, 33279, 62114,
                                                                       16303, 16393, 34929,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 64214, 0, 3,
                                                                       62114, 33379, 62264,
                                                                       16393, 16483, 35079,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 64439, 0, 3,
                                                                       62264, 33479, 62414,
                                                                       16483, 16573, 35229,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 64664, 0, 3,
                                                                       62414, 33579, 62564,
                                                                       16573, 16663, 35379,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 64889, 0, 3,
                                                                       62564, 33679, 62714,
                                                                       16663, 16753, 35529,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65114, 0, 3,
                                                                       62864, 34179, 63089,
                                                                       16933, 17059, 36099,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65429, 0, 3,
                                                                       63089, 34329, 63314,
                                                                       17059, 17185, 36309,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65744, 0, 3,
                                                                       63314, 34479, 63539,
                                                                       17185, 17311, 36519,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 66059, 0, 3,
                                                                       63539, 34629, 63764,
                                                                       17311, 17437, 36729,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 66374, 0, 3,
                                                                       63764, 34779, 63989,
                                                                       17437, 17563, 36939,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 66689, 0, 3,
                                                                       63989, 34929, 64214,
                                                                       17563, 17689, 37149,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 67004, 0, 3,
                                                                       64214, 35079, 64439,
                                                                       17689, 17815, 37359,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 67319, 0, 3,
                                                                       64439, 35229, 64664,
                                                                       17815, 17941, 37569,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 67634, 0, 3,
                                                                       64664, 35379, 64889,
                                                                       17941, 18067, 37779,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 67949, 0, 3,
                                                                       65114, 36099, 65429,
                                                                       18319, 18487, 38549,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 68369, 0, 3,
                                                                       65429, 36309, 65744,
                                                                       18487, 18655, 38829,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 68789, 0, 3,
                                                                       65744, 36519, 66059,
                                                                       18655, 18823, 39109,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 69209, 0, 3,
                                                                       66059, 36729, 66374,
                                                                       18823, 18991, 39389,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 69629, 0, 3,
                                                                       66374, 36939, 66689,
                                                                       18991, 19159, 39669,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 70049, 0, 3,
                                                                       66689, 37149, 67004,
                                                                       19159, 19327, 39949,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 70469, 0, 3,
                                                                       67004, 37359, 67319,
                                                                       19327, 19495, 40229,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 70889, 0, 3,
                                                                       67319, 37569, 67634,
                                                                       19495, 19663, 40509,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 71309, 0, 3,
                                                                       67949, 38549, 68369,
                                                                       19999, 20215, 41509,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 71849, 0, 3,
                                                                       68369, 38829, 68789,
                                                                       20215, 20431, 41869,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 72389, 0, 3,
                                                                       68789, 39109, 69209,
                                                                       20431, 20647, 42229,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 72929, 0, 3,
                                                                       69209, 39389, 69629,
                                                                       20647, 20863, 42589,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 73469, 0, 3,
                                                                       69629, 39669, 70049,
                                                                       20863, 21079, 42949,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 74009, 0, 3,
                                                                       70049, 39949, 70469,
                                                                       21079, 21295, 43309,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 74549, 0, 3,
                                                                       70469, 40229, 70889,
                                                                       21295, 21511, 43669,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 75089, 0, 3,
                                                                       71309, 41509, 71849,
                                                                       21943, 22213, 44929,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 75764, 0, 3,
                                                                       71849, 41869, 72389,
                                                                       22213, 22483, 45379,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 76439, 0, 3,
                                                                       72389, 42229, 72929,
                                                                       22483, 22753, 45829,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 77114, 0, 3,
                                                                       72929, 42589, 73469,
                                                                       22753, 23023, 46279,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 77789, 0, 3,
                                                                       73469, 42949, 74009,
                                                                       23023, 23293, 46729,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 78464, 0, 3,
                                                                       74009, 43309, 74549,
                                                                       23293, 23563, 47179,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 79139, 0, 3,
                                                                       75089, 44929, 75764,
                                                                       24103, 24433, 48729,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 79964, 0, 3,
                                                                       75764, 45379, 76439,
                                                                       24433, 24763, 49279,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 80789, 0, 3,
                                                                       76439, 45829, 77114,
                                                                       24763, 25093, 49829,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 81614, 0, 3,
                                                                       77114, 46279, 77789,
                                                                       25093, 25423, 50379,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 82439, 0, 3,
                                                                       77789, 46729, 78464,
                                                                       25423, 25753, 50929,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 83264, 0, 3,
                                                                       79139, 48729, 79964,
                                                                       26413, 26809, 52799,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 84254, 0, 3,
                                                                       79964, 49279, 80789,
                                                                       26809, 27205, 53459,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 85244, 0, 3,
                                                                       80789, 49829, 81614,
                                                                       27205, 27601, 54119,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 86234, 0, 3,
                                                                       81614, 50379, 82439,
                                                                       27601, 27997, 54779,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 87224, 0, 3,
                                                                       83264, 52799, 84254,
                                                                       28789, 29257, 56999,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 88394, 0, 3,
                                                                       84254, 53459, 85244,
                                                                       29257, 29725, 57779,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 89564, 0, 3,
                                                                       85244, 54119, 86234,
                                                                       29725, 30193, 58559,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90734, 3, 31129,
                                                                       31139, 59339, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90755, 3, 31139,
                                                                       31149, 59354, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90776, 3, 31149,
                                                                       31159, 59369, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90797, 3, 31159,
                                                                       31169, 59384, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90818, 3, 31169,
                                                                       31179, 59399, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90839, 3, 31179,
                                                                       31189, 59414, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90860, 3, 31189,
                                                                       31199, 59429, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90881, 3, 31199,
                                                                       31209, 59444, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90902, 3, 31209,
                                                                       31219, 59459, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90923, 3, 31219,
                                                                       31229, 59474, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90944, 3, 31229,
                                                                       31239, 59489, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90965, 3, 31239,
                                                                       31249, 59504, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90986, 3, 31249,
                                                                       31259, 59519, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 91007, 3, 31259,
                                                                       31269, 59534, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91028, 0, 3,
                                                                       90734, 59339, 90755,
                                                                       31289, 31319, 59549,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91091, 0, 3,
                                                                       90755, 59354, 90776,
                                                                       31319, 31349, 59594,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91154, 0, 3,
                                                                       90776, 59369, 90797,
                                                                       31349, 31379, 59639,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91217, 0, 3,
                                                                       90797, 59384, 90818,
                                                                       31379, 31409, 59684,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91280, 0, 3,
                                                                       90818, 59399, 90839,
                                                                       31409, 31439, 59729,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91343, 0, 3,
                                                                       90839, 59414, 90860,
                                                                       31439, 31469, 59774,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91406, 0, 3,
                                                                       90860, 59429, 90881,
                                                                       31469, 31499, 59819,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91469, 0, 3,
                                                                       90881, 59444, 90902,
                                                                       31499, 31529, 59864,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91532, 0, 3,
                                                                       90902, 59459, 90923,
                                                                       31529, 31559, 59909,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91595, 0, 3,
                                                                       90923, 59474, 90944,
                                                                       31559, 31589, 59954,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91658, 0, 3,
                                                                       90944, 59489, 90965,
                                                                       31589, 31619, 59999,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91721, 0, 3,
                                                                       90965, 59504, 90986,
                                                                       31619, 31649, 60044,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91784, 0, 3,
                                                                       90986, 59519, 91007,
                                                                       31649, 31679, 60089,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91847, 0, 3,
                                                                       91028, 59549, 91091,
                                                                       31739, 31799, 60134,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91973, 0, 3,
                                                                       91091, 59594, 91154,
                                                                       31799, 31859, 60224,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92099, 0, 3,
                                                                       91154, 59639, 91217,
                                                                       31859, 31919, 60314,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92225, 0, 3,
                                                                       91217, 59684, 91280,
                                                                       31919, 31979, 60404,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92351, 0, 3,
                                                                       91280, 59729, 91343,
                                                                       31979, 32039, 60494,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92477, 0, 3,
                                                                       91343, 59774, 91406,
                                                                       32039, 32099, 60584,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92603, 0, 3,
                                                                       91406, 59819, 91469,
                                                                       32099, 32159, 60674,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92729, 0, 3,
                                                                       91469, 59864, 91532,
                                                                       32159, 32219, 60764,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92855, 0, 3,
                                                                       91532, 59909, 91595,
                                                                       32219, 32279, 60854,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92981, 0, 3,
                                                                       91595, 59954, 91658,
                                                                       32279, 32339, 60944,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 93107, 0, 3,
                                                                       91658, 59999, 91721,
                                                                       32339, 32399, 61034,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 93233, 0, 3,
                                                                       91721, 60044, 91784,
                                                                       32399, 32459, 61124,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93359, 0, 3,
                                                                       91847, 60134, 91973,
                                                                       32579, 32679, 61214,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93569, 0, 3,
                                                                       91973, 60224, 92099,
                                                                       32679, 32779, 61364,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93779, 0, 3,
                                                                       92099, 60314, 92225,
                                                                       32779, 32879, 61514,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93989, 0, 3,
                                                                       92225, 60404, 92351,
                                                                       32879, 32979, 61664,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94199, 0, 3,
                                                                       92351, 60494, 92477,
                                                                       32979, 33079, 61814,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94409, 0, 3,
                                                                       92477, 60584, 92603,
                                                                       33079, 33179, 61964,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94619, 0, 3,
                                                                       92603, 60674, 92729,
                                                                       33179, 33279, 62114,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94829, 0, 3,
                                                                       92729, 60764, 92855,
                                                                       33279, 33379, 62264,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 95039, 0, 3,
                                                                       92855, 60854, 92981,
                                                                       33379, 33479, 62414,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 95249, 0, 3,
                                                                       92981, 60944, 93107,
                                                                       33479, 33579, 62564,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 95459, 0, 3,
                                                                       93107, 61034, 93233,
                                                                       33579, 33679, 62714,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 95669, 0, 3,
                                                                       93359, 61214, 93569,
                                                                       33879, 34029, 62864,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 95984, 0, 3,
                                                                       93569, 61364, 93779,
                                                                       34029, 34179, 63089,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 96299, 0, 3,
                                                                       93779, 61514, 93989,
                                                                       34179, 34329, 63314,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 96614, 0, 3,
                                                                       93989, 61664, 94199,
                                                                       34329, 34479, 63539,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 96929, 0, 3,
                                                                       94199, 61814, 94409,
                                                                       34479, 34629, 63764,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 97244, 0, 3,
                                                                       94409, 61964, 94619,
                                                                       34629, 34779, 63989,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 97559, 0, 3,
                                                                       94619, 62114, 94829,
                                                                       34779, 34929, 64214,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 97874, 0, 3,
                                                                       94829, 62264, 95039,
                                                                       34929, 35079, 64439,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 98189, 0, 3,
                                                                       95039, 62414, 95249,
                                                                       35079, 35229, 64664,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 98504, 0, 3,
                                                                       95249, 62564, 95459,
                                                                       35229, 35379, 64889,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 98819, 0, 3,
                                                                       95669, 62864, 95984,
                                                                       35679, 35889, 65114,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 99260, 0, 3,
                                                                       95984, 63089, 96299,
                                                                       35889, 36099, 65429,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 99701, 0, 3,
                                                                       96299, 63314, 96614,
                                                                       36099, 36309, 65744,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 100142, 0, 3,
                                                                       96614, 63539, 96929,
                                                                       36309, 36519, 66059,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 100583, 0, 3,
                                                                       96929, 63764, 97244,
                                                                       36519, 36729, 66374,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 101024, 0, 3,
                                                                       97244, 63989, 97559,
                                                                       36729, 36939, 66689,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 101465, 0, 3,
                                                                       97559, 64214, 97874,
                                                                       36939, 37149, 67004,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 101906, 0, 3,
                                                                       97874, 64439, 98189,
                                                                       37149, 37359, 67319,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 102347, 0, 3,
                                                                       98189, 64664, 98504,
                                                                       37359, 37569, 67634,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 102788, 0, 3,
                                                                       98819, 65114, 99260,
                                                                       37989, 38269, 67949,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 103376, 0, 3,
                                                                       99260, 65429, 99701,
                                                                       38269, 38549, 68369,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 103964, 0, 3,
                                                                       99701, 65744, 100142,
                                                                       38549, 38829, 68789,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 104552, 0, 3,
                                                                       100142, 66059, 100583,
                                                                       38829, 39109, 69209,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 105140, 0, 3,
                                                                       100583, 66374, 101024,
                                                                       39109, 39389, 69629,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 105728, 0, 3,
                                                                       101024, 66689, 101465,
                                                                       39389, 39669, 70049,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 106316, 0, 3,
                                                                       101465, 67004, 101906,
                                                                       39669, 39949, 70469,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 106904, 0, 3,
                                                                       101906, 67319, 102347,
                                                                       39949, 40229, 70889,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 107492, 0, 3,
                                                                       102788, 67949, 103376,
                                                                       40789, 41149, 71309,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 108248, 0, 3,
                                                                       103376, 68369, 103964,
                                                                       41149, 41509, 71849,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 109004, 0, 3,
                                                                       103964, 68789, 104552,
                                                                       41509, 41869, 72389,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 109760, 0, 3,
                                                                       104552, 69209, 105140,
                                                                       41869, 42229, 72929,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 110516, 0, 3,
                                                                       105140, 69629, 105728,
                                                                       42229, 42589, 73469,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 111272, 0, 3,
                                                                       105728, 70049, 106316,
                                                                       42589, 42949, 74009,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 112028, 0, 3,
                                                                       106316, 70469, 106904,
                                                                       42949, 43309, 74549,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 112784, 0, 3,
                                                                       107492, 71309, 108248,
                                                                       44029, 44479, 75089,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 113729, 0, 3,
                                                                       108248, 71849, 109004,
                                                                       44479, 44929, 75764,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 114674, 0, 3,
                                                                       109004, 72389, 109760,
                                                                       44929, 45379, 76439,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 115619, 0, 3,
                                                                       109760, 72929, 110516,
                                                                       45379, 45829, 77114,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 116564, 0, 3,
                                                                       110516, 73469, 111272,
                                                                       45829, 46279, 77789,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 117509, 0, 3,
                                                                       111272, 74009, 112028,
                                                                       46279, 46729, 78464,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 118454, 0, 3,
                                                                       112784, 75089, 113729,
                                                                       47629, 48179, 79139,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 119609, 0, 3,
                                                                       113729, 75764, 114674,
                                                                       48179, 48729, 79964,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 120764, 0, 3,
                                                                       114674, 76439, 115619,
                                                                       48729, 49279, 80789,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 121919, 0, 3,
                                                                       115619, 77114, 116564,
                                                                       49279, 49829, 81614,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 123074, 0, 3,
                                                                       116564, 77789, 117509,
                                                                       49829, 50379, 82439,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 124229, 0, 3,
                                                                       118454, 79139, 119609,
                                                                       51479, 52139, 83264,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 125615, 0, 3,
                                                                       119609, 79964, 120764,
                                                                       52139, 52799, 84254,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 127001, 0, 3,
                                                                       120764, 80789, 121919,
                                                                       52799, 53459, 85244,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 128387, 0, 3,
                                                                       121919, 81614, 123074,
                                                                       53459, 54119, 86234,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 129773, 0, 3,
                                                                       124229, 83264, 125615,
                                                                       55439, 56219, 87224,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 131411, 0, 3,
                                                                       125615, 84254, 127001,
                                                                       56219, 56999, 88394,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 133049, 0, 3,
                                                                       127001, 85244, 128387,
                                                                       56999, 57779, 89564,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134687, 3, 59339,
                                                                       59354, 90776, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134715, 3, 59354,
                                                                       59369, 90797, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134743, 3, 59369,
                                                                       59384, 90818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134771, 3, 59384,
                                                                       59399, 90839, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134799, 3, 59399,
                                                                       59414, 90860, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134827, 3, 59414,
                                                                       59429, 90881, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134855, 3, 59429,
                                                                       59444, 90902, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134883, 3, 59444,
                                                                       59459, 90923, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134911, 3, 59459,
                                                                       59474, 90944, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134939, 3, 59474,
                                                                       59489, 90965, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134967, 3, 59489,
                                                                       59504, 90986, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134995, 3, 59504,
                                                                       59519, 91007, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135023, 0, 3,
                                                                       134687, 90776, 134715,
                                                                       59549, 59594, 91154,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135107, 0, 3,
                                                                       134715, 90797, 134743,
                                                                       59594, 59639, 91217,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135191, 0, 3,
                                                                       134743, 90818, 134771,
                                                                       59639, 59684, 91280,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135275, 0, 3,
                                                                       134771, 90839, 134799,
                                                                       59684, 59729, 91343,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135359, 0, 3,
                                                                       134799, 90860, 134827,
                                                                       59729, 59774, 91406,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135443, 0, 3,
                                                                       134827, 90881, 134855,
                                                                       59774, 59819, 91469,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135527, 0, 3,
                                                                       134855, 90902, 134883,
                                                                       59819, 59864, 91532,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135611, 0, 3,
                                                                       134883, 90923, 134911,
                                                                       59864, 59909, 91595,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135695, 0, 3,
                                                                       134911, 90944, 134939,
                                                                       59909, 59954, 91658,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135779, 0, 3,
                                                                       134939, 90965, 134967,
                                                                       59954, 59999, 91721,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135863, 0, 3,
                                                                       134967, 90986, 134995,
                                                                       59999, 60044, 91784,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 135947, 0, 3,
                                                                       135023, 91154, 135107,
                                                                       60134, 60224, 92099,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136115, 0, 3,
                                                                       135107, 91217, 135191,
                                                                       60224, 60314, 92225,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136283, 0, 3,
                                                                       135191, 91280, 135275,
                                                                       60314, 60404, 92351,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136451, 0, 3,
                                                                       135275, 91343, 135359,
                                                                       60404, 60494, 92477,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136619, 0, 3,
                                                                       135359, 91406, 135443,
                                                                       60494, 60584, 92603,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136787, 0, 3,
                                                                       135443, 91469, 135527,
                                                                       60584, 60674, 92729,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136955, 0, 3,
                                                                       135527, 91532, 135611,
                                                                       60674, 60764, 92855,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 137123, 0, 3,
                                                                       135611, 91595, 135695,
                                                                       60764, 60854, 92981,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 137291, 0, 3,
                                                                       135695, 91658, 135779,
                                                                       60854, 60944, 93107,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 137459, 0, 3,
                                                                       135779, 91721, 135863,
                                                                       60944, 61034, 93233,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 137627, 0, 3,
                                                                       135947, 92099, 136115,
                                                                       61214, 61364, 93779,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 137907, 0, 3,
                                                                       136115, 92225, 136283,
                                                                       61364, 61514, 93989,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 138187, 0, 3,
                                                                       136283, 92351, 136451,
                                                                       61514, 61664, 94199,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 138467, 0, 3,
                                                                       136451, 92477, 136619,
                                                                       61664, 61814, 94409,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 138747, 0, 3,
                                                                       136619, 92603, 136787,
                                                                       61814, 61964, 94619,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 139027, 0, 3,
                                                                       136787, 92729, 136955,
                                                                       61964, 62114, 94829,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 139307, 0, 3,
                                                                       136955, 92855, 137123,
                                                                       62114, 62264, 95039,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 139587, 0, 3,
                                                                       137123, 92981, 137291,
                                                                       62264, 62414, 95249,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 139867, 0, 3,
                                                                       137291, 93107, 137459,
                                                                       62414, 62564, 95459,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 140147, 0, 3,
                                                                       137627, 93779, 137907,
                                                                       62864, 63089, 96299,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 140567, 0, 3,
                                                                       137907, 93989, 138187,
                                                                       63089, 63314, 96614,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 140987, 0, 3,
                                                                       138187, 94199, 138467,
                                                                       63314, 63539, 96929,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 141407, 0, 3,
                                                                       138467, 94409, 138747,
                                                                       63539, 63764, 97244,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 141827, 0, 3,
                                                                       138747, 94619, 139027,
                                                                       63764, 63989, 97559,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 142247, 0, 3,
                                                                       139027, 94829, 139307,
                                                                       63989, 64214, 97874,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 142667, 0, 3,
                                                                       139307, 95039, 139587,
                                                                       64214, 64439, 98189,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 143087, 0, 3,
                                                                       139587, 95249, 139867,
                                                                       64439, 64664, 98504,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 143507, 0, 3,
                                                                       140147, 96299, 140567,
                                                                       65114, 65429, 99701,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 144095, 0, 3,
                                                                       140567, 96614, 140987,
                                                                       65429, 65744, 100142,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 144683, 0, 3,
                                                                       140987, 96929, 141407,
                                                                       65744, 66059, 100583,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 145271, 0, 3,
                                                                       141407, 97244, 141827,
                                                                       66059, 66374, 101024,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 145859, 0, 3,
                                                                       141827, 97559, 142247,
                                                                       66374, 66689, 101465,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 146447, 0, 3,
                                                                       142247, 97874, 142667,
                                                                       66689, 67004, 101906,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 147035, 0, 3,
                                                                       142667, 98189, 143087,
                                                                       67004, 67319, 102347,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 147623, 0, 3,
                                                                       143507, 99701, 144095,
                                                                       67949, 68369, 103964,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 148407, 0, 3,
                                                                       144095, 100142, 144683,
                                                                       68369, 68789, 104552,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 149191, 0, 3,
                                                                       144683, 100583, 145271,
                                                                       68789, 69209, 105140,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 149975, 0, 3,
                                                                       145271, 101024, 145859,
                                                                       69209, 69629, 105728,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 150759, 0, 3,
                                                                       145859, 101465, 146447,
                                                                       69629, 70049, 106316,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 151543, 0, 3,
                                                                       146447, 101906, 147035,
                                                                       70049, 70469, 106904,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 152327, 0, 3,
                                                                       147623, 103964, 148407,
                                                                       71309, 71849, 109004,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 153335, 0, 3,
                                                                       148407, 104552, 149191,
                                                                       71849, 72389, 109760,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 154343, 0, 3,
                                                                       149191, 105140, 149975,
                                                                       72389, 72929, 110516,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 155351, 0, 3,
                                                                       149975, 105728, 150759,
                                                                       72929, 73469, 111272,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 156359, 0, 3,
                                                                       150759, 106316, 151543,
                                                                       73469, 74009, 112028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 157367, 0, 3,
                                                                       152327, 109004, 153335,
                                                                       75089, 75764, 114674,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 158627, 0, 3,
                                                                       153335, 109760, 154343,
                                                                       75764, 76439, 115619,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 159887, 0, 3,
                                                                       154343, 110516, 155351,
                                                                       76439, 77114, 116564,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 161147, 0, 3,
                                                                       155351, 111272, 156359,
                                                                       77114, 77789, 117509,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 162407, 0, 3,
                                                                       157367, 114674, 158627,
                                                                       79139, 79964, 120764,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 163947, 0, 3,
                                                                       158627, 115619, 159887,
                                                                       79964, 80789, 121919,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 165487, 0, 3,
                                                                       159887, 116564, 161147,
                                                                       80789, 81614, 123074,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 167027, 0, 3,
                                                                       162407, 120764, 163947,
                                                                       83264, 84254, 127001,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 168875, 0, 3,
                                                                       163947, 121919, 165487,
                                                                       84254, 85244, 128387,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 170723, 0, 3,
                                                                       167027, 127001, 168875,
                                                                       87224, 88394, 133049,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172907, 3, 90734,
                                                                       90755, 134687, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172943, 3, 90755,
                                                                       90776, 134715, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172979, 3, 90776,
                                                                       90797, 134743, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173015, 3, 90797,
                                                                       90818, 134771, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173051, 3, 90818,
                                                                       90839, 134799, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173087, 3, 90839,
                                                                       90860, 134827, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173123, 3, 90860,
                                                                       90881, 134855, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173159, 3, 90881,
                                                                       90902, 134883, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173195, 3, 90902,
                                                                       90923, 134911, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173231, 3, 90923,
                                                                       90944, 134939, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173267, 3, 90944,
                                                                       90965, 134967, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173303, 3, 90965,
                                                                       90986, 134995, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173339, 0, 3,
                                                                       172907, 134687, 172943,
                                                                       91028, 91091, 135023,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173447, 0, 3,
                                                                       172943, 134715, 172979,
                                                                       91091, 91154, 135107,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173555, 0, 3,
                                                                       172979, 134743, 173015,
                                                                       91154, 91217, 135191,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173663, 0, 3,
                                                                       173015, 134771, 173051,
                                                                       91217, 91280, 135275,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173771, 0, 3,
                                                                       173051, 134799, 173087,
                                                                       91280, 91343, 135359,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173879, 0, 3,
                                                                       173087, 134827, 173123,
                                                                       91343, 91406, 135443,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173987, 0, 3,
                                                                       173123, 134855, 173159,
                                                                       91406, 91469, 135527,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 174095, 0, 3,
                                                                       173159, 134883, 173195,
                                                                       91469, 91532, 135611,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 174203, 0, 3,
                                                                       173195, 134911, 173231,
                                                                       91532, 91595, 135695,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 174311, 0, 3,
                                                                       173231, 134939, 173267,
                                                                       91595, 91658, 135779,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 174419, 0, 3,
                                                                       173267, 134967, 173303,
                                                                       91658, 91721, 135863,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 174527, 0, 3,
                                                                       173339, 135023, 173447,
                                                                       91847, 91973, 135947,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 174743, 0, 3,
                                                                       173447, 135107, 173555,
                                                                       91973, 92099, 136115,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 174959, 0, 3,
                                                                       173555, 135191, 173663,
                                                                       92099, 92225, 136283,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 175175, 0, 3,
                                                                       173663, 135275, 173771,
                                                                       92225, 92351, 136451,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 175391, 0, 3,
                                                                       173771, 135359, 173879,
                                                                       92351, 92477, 136619,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 175607, 0, 3,
                                                                       173879, 135443, 173987,
                                                                       92477, 92603, 136787,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 175823, 0, 3,
                                                                       173987, 135527, 174095,
                                                                       92603, 92729, 136955,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 176039, 0, 3,
                                                                       174095, 135611, 174203,
                                                                       92729, 92855, 137123,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 176255, 0, 3,
                                                                       174203, 135695, 174311,
                                                                       92855, 92981, 137291,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 176471, 0, 3,
                                                                       174311, 135779, 174419,
                                                                       92981, 93107, 137459,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 176687, 0, 3,
                                                                       174527, 135947, 174743,
                                                                       93359, 93569, 137627,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 177047, 0, 3,
                                                                       174743, 136115, 174959,
                                                                       93569, 93779, 137907,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 177407, 0, 3,
                                                                       174959, 136283, 175175,
                                                                       93779, 93989, 138187,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 177767, 0, 3,
                                                                       175175, 136451, 175391,
                                                                       93989, 94199, 138467,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 178127, 0, 3,
                                                                       175391, 136619, 175607,
                                                                       94199, 94409, 138747,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 178487, 0, 3,
                                                                       175607, 136787, 175823,
                                                                       94409, 94619, 139027,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 178847, 0, 3,
                                                                       175823, 136955, 176039,
                                                                       94619, 94829, 139307,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 179207, 0, 3,
                                                                       176039, 137123, 176255,
                                                                       94829, 95039, 139587,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 179567, 0, 3,
                                                                       176255, 137291, 176471,
                                                                       95039, 95249, 139867,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 179927, 0, 3,
                                                                       176687, 137627, 177047,
                                                                       95669, 95984, 140147,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 180467, 0, 3,
                                                                       177047, 137907, 177407,
                                                                       95984, 96299, 140567,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 181007, 0, 3,
                                                                       177407, 138187, 177767,
                                                                       96299, 96614, 140987,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 181547, 0, 3,
                                                                       177767, 138467, 178127,
                                                                       96614, 96929, 141407,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 182087, 0, 3,
                                                                       178127, 138747, 178487,
                                                                       96929, 97244, 141827,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 182627, 0, 3,
                                                                       178487, 139027, 178847,
                                                                       97244, 97559, 142247,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 183167, 0, 3,
                                                                       178847, 139307, 179207,
                                                                       97559, 97874, 142667,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 183707, 0, 3,
                                                                       179207, 139587, 179567,
                                                                       97874, 98189, 143087,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 184247, 0, 3,
                                                                       179927, 140147, 180467,
                                                                       98819, 99260, 143507,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 185003, 0, 3,
                                                                       180467, 140567, 181007,
                                                                       99260, 99701, 144095,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 185759, 0, 3,
                                                                       181007, 140987, 181547,
                                                                       99701, 100142, 144683,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 186515, 0, 3,
                                                                       181547, 141407, 182087,
                                                                       100142, 100583, 145271,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 187271, 0, 3,
                                                                       182087, 141827, 182627,
                                                                       100583, 101024, 145859,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 188027, 0, 3,
                                                                       182627, 142247, 183167,
                                                                       101024, 101465, 146447,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 188783, 0, 3,
                                                                       183167, 142667, 183707,
                                                                       101465, 101906, 147035,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 189539, 0, 3,
                                                                       184247, 143507, 185003,
                                                                       102788, 103376, 147623,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 190547, 0, 3,
                                                                       185003, 144095, 185759,
                                                                       103376, 103964, 148407,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 191555, 0, 3,
                                                                       185759, 144683, 186515,
                                                                       103964, 104552, 149191,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 192563, 0, 3,
                                                                       186515, 145271, 187271,
                                                                       104552, 105140, 149975,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 193571, 0, 3,
                                                                       187271, 145859, 188027,
                                                                       105140, 105728, 150759,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 194579, 0, 3,
                                                                       188027, 146447, 188783,
                                                                       105728, 106316, 151543,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 195587, 0, 3,
                                                                       189539, 147623, 190547,
                                                                       107492, 108248, 152327,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 196883, 0, 3,
                                                                       190547, 148407, 191555,
                                                                       108248, 109004, 153335,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 198179, 0, 3,
                                                                       191555, 149191, 192563,
                                                                       109004, 109760, 154343,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 199475, 0, 3,
                                                                       192563, 149975, 193571,
                                                                       109760, 110516, 155351,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 200771, 0, 3,
                                                                       193571, 150759, 194579,
                                                                       110516, 111272, 156359,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 202067, 0, 3,
                                                                       195587, 152327, 196883,
                                                                       112784, 113729, 157367,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 203687, 0, 3,
                                                                       196883, 153335, 198179,
                                                                       113729, 114674, 158627,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 205307, 0, 3,
                                                                       198179, 154343, 199475,
                                                                       114674, 115619, 159887,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 206927, 0, 3,
                                                                       199475, 155351, 200771,
                                                                       115619, 116564, 161147,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 208547, 0, 3,
                                                                       202067, 157367, 203687,
                                                                       118454, 119609, 162407,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 210527, 0, 3,
                                                                       203687, 158627, 205307,
                                                                       119609, 120764, 163947,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 212507, 0, 3,
                                                                       205307, 159887, 206927,
                                                                       120764, 121919, 165487,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 214487, 0, 3,
                                                                       208547, 162407, 210527,
                                                                       124229, 125615, 167027,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 216863, 0, 3,
                                                                       210527, 163947, 212507,
                                                                       125615, 127001, 168875,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 219239, 0, 3,
                                                                       214487, 167027, 216863,
                                                                       129773, 131411, 170723,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 222047, 189539, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 223475, 195587, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 225311, 202067, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 227606, 208547, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 230411, 214487, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 233777, 219239, 2808, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 223055, 222047, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 224771, 223475, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 226931, 225311, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 229586, 227606, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 232787, 230411, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 236585, 233777, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 237755, 223055, 224771, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 239015, 224771, 226931, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 240635, 226931, 229586, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 242660, 229586, 232787, 15,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 245135, 232787, 236585, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 248105, 237755, 239015, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 250625, 239015, 240635, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 253865, 240635, 242660, 15,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 257915, 242660, 245135, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 262865, 248105, 250625, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 267065, 250625, 253865, 15,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 272465, 253865, 257915, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 279215, 262865, 267065, 15,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 285515, 267065, 272465, 15,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 293615, 279215, 285515, 15,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 302435, 293615, 28, 15, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 302435, 165, nmax);
    }

    for (size_t m = 0; m < 2145; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
