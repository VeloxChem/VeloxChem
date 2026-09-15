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

    const auto nmax = simdfunc::prepare_buffer(buffer, 307056, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 307056, 222048, 14538, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 7, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16, 17, 18}, ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 8, 9,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 9, 10,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 10, 11,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 11, 12,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 101, 0, 3, 12, 13,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 107, 0, 3, 13, 14,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 14, 15,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 119, 0, 3, 15, 16,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 125, 0, 3, 16, 17,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 131, 0, 3, 17, 18,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 137, 0, 3, 18, 19,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 19, 20,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 149, 0, 3, 20, 21,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 155, 0, 3, 21, 22,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 161, 0, 3, 22, 23,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 167, 0, 3, 23, 24,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 26, 29,
                                                                       77, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 29, 32,
                                                                       83, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 32, 35,
                                                                       89, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 35, 38,
                                                                       95, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 38, 41,
                                                                       101, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 41, 44,
                                                                       107, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 44, 47,
                                                                       113, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 47, 50,
                                                                       119, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 50, 53,
                                                                       125, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 263, 0, 3, 53, 56,
                                                                       131, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 56, 59,
                                                                       137, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 283, 0, 3, 59, 62,
                                                                       143, 149, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 293, 0, 3, 62, 65,
                                                                       149, 155, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 303, 0, 3, 65, 68,
                                                                       155, 161, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 313, 0, 3, 68, 71,
                                                                       161, 167, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 323, 0, 3, 77, 83,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 83, 89,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 353, 0, 3, 89, 95,
                                                                       193, 203, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 95,
                                                                       101, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 383, 0, 3, 101,
                                                                       107, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 107,
                                                                       113, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 413, 0, 3, 113,
                                                                       119, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 119,
                                                                       125, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 443, 0, 3, 125,
                                                                       131, 253, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 131,
                                                                       137, 263, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 473, 0, 3, 137,
                                                                       143, 273, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 143,
                                                                       149, 283, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 503, 0, 3, 149,
                                                                       155, 293, 303, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 155,
                                                                       161, 303, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 533, 0, 3, 173,
                                                                       183, 323, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 554, 0, 3, 183,
                                                                       193, 338, 353, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 575, 0, 3, 193,
                                                                       203, 353, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 596, 0, 3, 203,
                                                                       213, 368, 383, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 617, 0, 3, 213,
                                                                       223, 383, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 223,
                                                                       233, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 659, 0, 3, 233,
                                                                       243, 413, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 680, 0, 3, 243,
                                                                       253, 428, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 701, 0, 3, 253,
                                                                       263, 443, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 722, 0, 3, 263,
                                                                       273, 458, 473, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 743, 0, 3, 273,
                                                                       283, 473, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 764, 0, 3, 283,
                                                                       293, 488, 503, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 785, 0, 3, 293,
                                                                       303, 503, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 806, 0, 3, 323,
                                                                       338, 533, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 834, 0, 3, 338,
                                                                       353, 554, 575, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 862, 0, 3, 353,
                                                                       368, 575, 596, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 890, 0, 3, 368,
                                                                       383, 596, 617, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 918, 0, 3, 383,
                                                                       398, 617, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 946, 0, 3, 398,
                                                                       413, 638, 659, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 974, 0, 3, 413,
                                                                       428, 659, 680, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 428,
                                                                       443, 680, 701, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 443,
                                                                       458, 701, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 458,
                                                                       473, 722, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1086, 0, 3, 473,
                                                                       488, 743, 764, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 488,
                                                                       503, 764, 785, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 533,
                                                                       554, 806, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1178, 0, 3, 554,
                                                                       575, 834, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1214, 0, 3, 575,
                                                                       596, 862, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1250, 0, 3, 596,
                                                                       617, 890, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1286, 0, 3, 617,
                                                                       638, 918, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1322, 0, 3, 638,
                                                                       659, 946, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1358, 0, 3, 659,
                                                                       680, 974, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1394, 0, 3, 680,
                                                                       701, 1002, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1430, 0, 3, 701,
                                                                       722, 1030, 1058, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1466, 0, 3, 722,
                                                                       743, 1058, 1086, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1502, 0, 3, 743,
                                                                       764, 1086, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1538, 0, 3, 806,
                                                                       834, 1142, 1178, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1583, 0, 3, 834,
                                                                       862, 1178, 1214, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1628, 0, 3, 862,
                                                                       890, 1214, 1250, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1673, 0, 3, 890,
                                                                       918, 1250, 1286, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1718, 0, 3, 918,
                                                                       946, 1286, 1322, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1763, 0, 3, 946,
                                                                       974, 1322, 1358, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1808, 0, 3, 974,
                                                                       1002, 1358, 1394, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1853, 0, 3, 1002,
                                                                       1030, 1394, 1430, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1898, 0, 3, 1030,
                                                                       1058, 1430, 1466, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1943, 0, 3, 1058,
                                                                       1086, 1466, 1502, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1142,
                                                                       1178, 1538, 1583, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2043, 0, 3, 1178,
                                                                       1214, 1583, 1628, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2098, 0, 3, 1214,
                                                                       1250, 1628, 1673, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1250,
                                                                       1286, 1673, 1718, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2208, 0, 3, 1286,
                                                                       1322, 1718, 1763, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2263, 0, 3, 1322,
                                                                       1358, 1763, 1808, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2318, 0, 3, 1358,
                                                                       1394, 1808, 1853, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2373, 0, 3, 1394,
                                                                       1430, 1853, 1898, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2428, 0, 3, 1430,
                                                                       1466, 1898, 1943, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2483, 0, 3, 1538,
                                                                       1583, 1988, 2043, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2549, 0, 3, 1583,
                                                                       1628, 2043, 2098, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2615, 0, 3, 1628,
                                                                       1673, 2098, 2153, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2681, 0, 3, 1673,
                                                                       1718, 2153, 2208, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2747, 0, 3, 1718,
                                                                       1763, 2208, 2263, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2813, 0, 3, 1763,
                                                                       1808, 2263, 2318, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2879, 0, 3, 1808,
                                                                       1853, 2318, 2373, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2945, 0, 3, 1853,
                                                                       1898, 2373, 2428, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3011, 0, 3, 1988,
                                                                       2043, 2483, 2549, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3089, 0, 3, 2043,
                                                                       2098, 2549, 2615, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3167, 0, 3, 2098,
                                                                       2153, 2615, 2681, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3245, 0, 3, 2153,
                                                                       2208, 2681, 2747, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3323, 0, 3, 2208,
                                                                       2263, 2747, 2813, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3401, 0, 3, 2263,
                                                                       2318, 2813, 2879, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3479, 0, 3, 2318,
                                                                       2373, 2879, 2945, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3557, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3560, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3563, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3566, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3569, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3572, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3575, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3578, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3581, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3584, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3587, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3590, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3593, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3596, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3599, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3602, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3605, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3608, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3611, 3, 8, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3620, 3, 9, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3629, 3, 10, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3638, 3, 11, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3647, 3, 12, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3656, 3, 13, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3665, 3, 14, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3674, 3, 15, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3683, 3, 16, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3692, 3, 17, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3701, 3, 18, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3710, 3, 19, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3719, 3, 20, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3728, 3, 21, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3737, 3, 22, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3746, 3, 23, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3755, 3, 24, 74,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3764, 3, 26, 77,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3782, 3, 29, 83,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3800, 3, 32, 89,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3818, 3, 35, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3836, 3, 38, 101,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3854, 3, 41, 107,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3872, 3, 44, 113,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3890, 3, 47, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3908, 3, 50, 125,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3926, 3, 53, 131,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3944, 3, 56, 137,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3962, 3, 59, 143,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3980, 3, 62, 149,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3998, 3, 65, 155,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4016, 3, 68, 161,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4034, 3, 71, 167,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4052, 3, 77, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4082, 3, 83, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4112, 3, 89, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4142, 3, 95, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4172, 3, 101, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4202, 3, 107, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4232, 3, 113, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4262, 3, 119, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4292, 3, 125, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4322, 3, 131, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4352, 3, 137, 273,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4382, 3, 143, 283,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4412, 3, 149, 293,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4442, 3, 155, 303,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4472, 3, 161, 313,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4502, 3, 173, 323,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4547, 3, 183, 338,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4592, 3, 193, 353,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4637, 3, 203, 368,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4682, 3, 213, 383,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4727, 3, 223, 398,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4772, 3, 233, 413,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4817, 3, 243, 428,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4862, 3, 253, 443,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4907, 3, 263, 458,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4952, 3, 273, 473,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4997, 3, 283, 488,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5042, 3, 293, 503,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5087, 3, 303, 518,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5132, 3, 323, 533,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5195, 3, 338, 554,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5258, 3, 353, 575,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5321, 3, 368, 596,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5384, 3, 383, 617,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5447, 3, 398, 638,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5510, 3, 413, 659,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5573, 3, 428, 680,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5636, 3, 443, 701,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5699, 3, 458, 722,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5762, 3, 473, 743,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5825, 3, 488, 764,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5888, 3, 503, 785,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5951, 3, 533, 806,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6035, 3, 554, 834,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6119, 3, 575, 862,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6203, 3, 596, 890,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6287, 3, 617, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6371, 3, 638, 946,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6455, 3, 659, 974,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6539, 3, 680,
                                                                       1002, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6623, 3, 701,
                                                                       1030, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6707, 3, 722,
                                                                       1058, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6791, 3, 743,
                                                                       1086, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6875, 3, 764,
                                                                       1114, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6959, 3, 806,
                                                                       1142, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7067, 3, 834,
                                                                       1178, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7175, 3, 862,
                                                                       1214, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7283, 3, 890,
                                                                       1250, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7391, 3, 918,
                                                                       1286, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7499, 3, 946,
                                                                       1322, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7607, 3, 974,
                                                                       1358, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7715, 3, 1002,
                                                                       1394, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7823, 3, 1030,
                                                                       1430, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7931, 3, 1058,
                                                                       1466, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8039, 3, 1086,
                                                                       1502, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8147, 3, 1142,
                                                                       1538, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8282, 3, 1178,
                                                                       1583, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8417, 3, 1214,
                                                                       1628, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8552, 3, 1250,
                                                                       1673, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8687, 3, 1286,
                                                                       1718, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8822, 3, 1322,
                                                                       1763, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8957, 3, 1358,
                                                                       1808, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9092, 3, 1394,
                                                                       1853, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9227, 3, 1430,
                                                                       1898, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9362, 3, 1466,
                                                                       1943, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9497, 3, 1538,
                                                                       1988, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9662, 3, 1583,
                                                                       2043, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9827, 3, 1628,
                                                                       2098, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9992, 3, 1673,
                                                                       2153, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10157, 3, 1718,
                                                                       2208, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10322, 3, 1763,
                                                                       2263, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10487, 3, 1808,
                                                                       2318, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10652, 3, 1853,
                                                                       2373, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10817, 3, 1898,
                                                                       2428, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10982, 3, 1988,
                                                                       2483, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11180, 3, 2043,
                                                                       2549, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11378, 3, 2098,
                                                                       2615, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11576, 3, 2153,
                                                                       2681, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11774, 3, 2208,
                                                                       2747, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11972, 3, 2263,
                                                                       2813, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12170, 3, 2318,
                                                                       2879, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12368, 3, 2373,
                                                                       2945, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12566, 3, 2483,
                                                                       3011, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12800, 3, 2549,
                                                                       3089, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13034, 3, 2615,
                                                                       3167, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13268, 3, 2681,
                                                                       3245, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13502, 3, 2747,
                                                                       3323, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13736, 3, 2813,
                                                                       3401, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13970, 3, 2879,
                                                                       3479, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14204, 3, 8, 9,
                                                                       3563, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14210, 3, 9, 10,
                                                                       3566, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14216, 3, 10, 11,
                                                                       3569, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14222, 3, 11, 12,
                                                                       3572, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14228, 3, 12, 13,
                                                                       3575, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14234, 3, 13, 14,
                                                                       3578, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14240, 3, 14, 15,
                                                                       3581, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14246, 3, 15, 16,
                                                                       3584, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14252, 3, 16, 17,
                                                                       3587, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14258, 3, 17, 18,
                                                                       3590, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14264, 3, 18, 19,
                                                                       3593, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14270, 3, 19, 20,
                                                                       3596, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14276, 3, 20, 21,
                                                                       3599, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14282, 3, 21, 22,
                                                                       3602, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14288, 3, 22, 23,
                                                                       3605, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14294, 3, 23, 24,
                                                                       3608, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14300, 0, 3,
                                                                       14204, 3563, 14210, 3629,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14318, 0, 3,
                                                                       14210, 3566, 14216, 3638,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14336, 0, 3,
                                                                       14216, 3569, 14222, 3647,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14354, 0, 3,
                                                                       14222, 3572, 14228, 3656,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14372, 0, 3,
                                                                       14228, 3575, 14234, 3665,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14390, 0, 3,
                                                                       14234, 3578, 14240, 3674,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14408, 0, 3,
                                                                       14240, 3581, 14246, 3683,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14426, 0, 3,
                                                                       14246, 3584, 14252, 3692,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14444, 0, 3,
                                                                       14252, 3587, 14258, 3701,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14462, 0, 3,
                                                                       14258, 3590, 14264, 3710,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14480, 0, 3,
                                                                       14264, 3593, 14270, 3719,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14498, 0, 3,
                                                                       14270, 3596, 14276, 3728,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14516, 0, 3,
                                                                       14276, 3599, 14282, 3737,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14534, 0, 3,
                                                                       14282, 3602, 14288, 3746,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14552, 0, 3,
                                                                       14288, 3605, 14294, 3755,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14570, 0, 3,
                                                                       14300, 3629, 14318, 77,
                                                                       83, 3800, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14606, 0, 3,
                                                                       14318, 3638, 14336, 83,
                                                                       89, 3818, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14642, 0, 3,
                                                                       14336, 3647, 14354, 89,
                                                                       95, 3836, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14678, 0, 3,
                                                                       14354, 3656, 14372, 95,
                                                                       101, 3854, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14714, 0, 3,
                                                                       14372, 3665, 14390, 101,
                                                                       107, 3872, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14750, 0, 3,
                                                                       14390, 3674, 14408, 107,
                                                                       113, 3890, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14786, 0, 3,
                                                                       14408, 3683, 14426, 113,
                                                                       119, 3908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14822, 0, 3,
                                                                       14426, 3692, 14444, 119,
                                                                       125, 3926, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14858, 0, 3,
                                                                       14444, 3701, 14462, 125,
                                                                       131, 3944, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14894, 0, 3,
                                                                       14462, 3710, 14480, 131,
                                                                       137, 3962, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14930, 0, 3,
                                                                       14480, 3719, 14498, 137,
                                                                       143, 3980, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14966, 0, 3,
                                                                       14498, 3728, 14516, 143,
                                                                       149, 3998, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15002, 0, 3,
                                                                       14516, 3737, 14534, 149,
                                                                       155, 4016, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15038, 0, 3,
                                                                       14534, 3746, 14552, 155,
                                                                       161, 4034, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15074, 0, 3,
                                                                       14570, 3800, 14606, 173,
                                                                       183, 4112, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15134, 0, 3,
                                                                       14606, 3818, 14642, 183,
                                                                       193, 4142, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15194, 0, 3,
                                                                       14642, 3836, 14678, 193,
                                                                       203, 4172, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15254, 0, 3,
                                                                       14678, 3854, 14714, 203,
                                                                       213, 4202, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15314, 0, 3,
                                                                       14714, 3872, 14750, 213,
                                                                       223, 4232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15374, 0, 3,
                                                                       14750, 3890, 14786, 223,
                                                                       233, 4262, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15434, 0, 3,
                                                                       14786, 3908, 14822, 233,
                                                                       243, 4292, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15494, 0, 3,
                                                                       14822, 3926, 14858, 243,
                                                                       253, 4322, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15554, 0, 3,
                                                                       14858, 3944, 14894, 253,
                                                                       263, 4352, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15614, 0, 3,
                                                                       14894, 3962, 14930, 263,
                                                                       273, 4382, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15674, 0, 3,
                                                                       14930, 3980, 14966, 273,
                                                                       283, 4412, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15734, 0, 3,
                                                                       14966, 3998, 15002, 283,
                                                                       293, 4442, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15794, 0, 3,
                                                                       15002, 4016, 15038, 293,
                                                                       303, 4472, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15854, 0, 3,
                                                                       15074, 4112, 15134, 323,
                                                                       338, 4592, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15944, 0, 3,
                                                                       15134, 4142, 15194, 338,
                                                                       353, 4637, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16034, 0, 3,
                                                                       15194, 4172, 15254, 353,
                                                                       368, 4682, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16124, 0, 3,
                                                                       15254, 4202, 15314, 368,
                                                                       383, 4727, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16214, 0, 3,
                                                                       15314, 4232, 15374, 383,
                                                                       398, 4772, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16304, 0, 3,
                                                                       15374, 4262, 15434, 398,
                                                                       413, 4817, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16394, 0, 3,
                                                                       15434, 4292, 15494, 413,
                                                                       428, 4862, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16484, 0, 3,
                                                                       15494, 4322, 15554, 428,
                                                                       443, 4907, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16574, 0, 3,
                                                                       15554, 4352, 15614, 443,
                                                                       458, 4952, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16664, 0, 3,
                                                                       15614, 4382, 15674, 458,
                                                                       473, 4997, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16754, 0, 3,
                                                                       15674, 4412, 15734, 473,
                                                                       488, 5042, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16844, 0, 3,
                                                                       15734, 4442, 15794, 488,
                                                                       503, 5087, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16934, 0, 3,
                                                                       15854, 4592, 15944, 533,
                                                                       554, 5258, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17060, 0, 3,
                                                                       15944, 4637, 16034, 554,
                                                                       575, 5321, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17186, 0, 3,
                                                                       16034, 4682, 16124, 575,
                                                                       596, 5384, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17312, 0, 3,
                                                                       16124, 4727, 16214, 596,
                                                                       617, 5447, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17438, 0, 3,
                                                                       16214, 4772, 16304, 617,
                                                                       638, 5510, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17564, 0, 3,
                                                                       16304, 4817, 16394, 638,
                                                                       659, 5573, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17690, 0, 3,
                                                                       16394, 4862, 16484, 659,
                                                                       680, 5636, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17816, 0, 3,
                                                                       16484, 4907, 16574, 680,
                                                                       701, 5699, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17942, 0, 3,
                                                                       16574, 4952, 16664, 701,
                                                                       722, 5762, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18068, 0, 3,
                                                                       16664, 4997, 16754, 722,
                                                                       743, 5825, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18194, 0, 3,
                                                                       16754, 5042, 16844, 743,
                                                                       764, 5888, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18320, 0, 3,
                                                                       16934, 5258, 17060, 806,
                                                                       834, 6119, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18488, 0, 3,
                                                                       17060, 5321, 17186, 834,
                                                                       862, 6203, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18656, 0, 3,
                                                                       17186, 5384, 17312, 862,
                                                                       890, 6287, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18824, 0, 3,
                                                                       17312, 5447, 17438, 890,
                                                                       918, 6371, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18992, 0, 3,
                                                                       17438, 5510, 17564, 918,
                                                                       946, 6455, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19160, 0, 3,
                                                                       17564, 5573, 17690, 946,
                                                                       974, 6539, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19328, 0, 3,
                                                                       17690, 5636, 17816, 974,
                                                                       1002, 6623, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19496, 0, 3,
                                                                       17816, 5699, 17942, 1002,
                                                                       1030, 6707, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19664, 0, 3,
                                                                       17942, 5762, 18068, 1030,
                                                                       1058, 6791, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19832, 0, 3,
                                                                       18068, 5825, 18194, 1058,
                                                                       1086, 6875, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20000, 0, 3,
                                                                       18320, 6119, 18488, 1142,
                                                                       1178, 7175, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20216, 0, 3,
                                                                       18488, 6203, 18656, 1178,
                                                                       1214, 7283, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20432, 0, 3,
                                                                       18656, 6287, 18824, 1214,
                                                                       1250, 7391, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20648, 0, 3,
                                                                       18824, 6371, 18992, 1250,
                                                                       1286, 7499, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20864, 0, 3,
                                                                       18992, 6455, 19160, 1286,
                                                                       1322, 7607, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21080, 0, 3,
                                                                       19160, 6539, 19328, 1322,
                                                                       1358, 7715, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21296, 0, 3,
                                                                       19328, 6623, 19496, 1358,
                                                                       1394, 7823, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21512, 0, 3,
                                                                       19496, 6707, 19664, 1394,
                                                                       1430, 7931, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21728, 0, 3,
                                                                       19664, 6791, 19832, 1430,
                                                                       1466, 8039, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21944, 0, 3,
                                                                       20000, 7175, 20216, 1538,
                                                                       1583, 8417, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22214, 0, 3,
                                                                       20216, 7283, 20432, 1583,
                                                                       1628, 8552, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22484, 0, 3,
                                                                       20432, 7391, 20648, 1628,
                                                                       1673, 8687, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22754, 0, 3,
                                                                       20648, 7499, 20864, 1673,
                                                                       1718, 8822, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23024, 0, 3,
                                                                       20864, 7607, 21080, 1718,
                                                                       1763, 8957, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23294, 0, 3,
                                                                       21080, 7715, 21296, 1763,
                                                                       1808, 9092, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23564, 0, 3,
                                                                       21296, 7823, 21512, 1808,
                                                                       1853, 9227, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23834, 0, 3,
                                                                       21512, 7931, 21728, 1853,
                                                                       1898, 9362, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24104, 0, 3,
                                                                       21944, 8417, 22214, 1988,
                                                                       2043, 9827, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24434, 0, 3,
                                                                       22214, 8552, 22484, 2043,
                                                                       2098, 9992, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24764, 0, 3,
                                                                       22484, 8687, 22754, 2098,
                                                                       2153, 10157, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25094, 0, 3,
                                                                       22754, 8822, 23024, 2153,
                                                                       2208, 10322, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25424, 0, 3,
                                                                       23024, 8957, 23294, 2208,
                                                                       2263, 10487, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25754, 0, 3,
                                                                       23294, 9092, 23564, 2263,
                                                                       2318, 10652, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 26084, 0, 3,
                                                                       23564, 9227, 23834, 2318,
                                                                       2373, 10817, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26414, 0, 3,
                                                                       24104, 9827, 24434, 2483,
                                                                       2549, 11378, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26810, 0, 3,
                                                                       24434, 9992, 24764, 2549,
                                                                       2615, 11576, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27206, 0, 3,
                                                                       24764, 10157, 25094, 2615,
                                                                       2681, 11774, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27602, 0, 3,
                                                                       25094, 10322, 25424, 2681,
                                                                       2747, 11972, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27998, 0, 3,
                                                                       25424, 10487, 25754, 2747,
                                                                       2813, 12170, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 28394, 0, 3,
                                                                       25754, 10652, 26084, 2813,
                                                                       2879, 12368, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 28790, 0, 3,
                                                                       26414, 11378, 26810, 3011,
                                                                       3089, 13034, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 29258, 0, 3,
                                                                       26810, 11576, 27206, 3089,
                                                                       3167, 13268, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 29726, 0, 3,
                                                                       27206, 11774, 27602, 3167,
                                                                       3245, 13502, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 30194, 0, 3,
                                                                       27602, 11972, 27998, 3245,
                                                                       3323, 13736, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 30662, 0, 3,
                                                                       27998, 12170, 28394, 3323,
                                                                       3401, 13970, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31130, 3, 3557,
                                                                       3560, 14204, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31140, 3, 3560,
                                                                       3563, 14210, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31150, 3, 3563,
                                                                       3566, 14216, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31160, 3, 3566,
                                                                       3569, 14222, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31170, 3, 3569,
                                                                       3572, 14228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31180, 3, 3572,
                                                                       3575, 14234, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31190, 3, 3575,
                                                                       3578, 14240, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31200, 3, 3578,
                                                                       3581, 14246, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31210, 3, 3581,
                                                                       3584, 14252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31220, 3, 3584,
                                                                       3587, 14258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31230, 3, 3587,
                                                                       3590, 14264, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31240, 3, 3590,
                                                                       3593, 14270, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31250, 3, 3593,
                                                                       3596, 14276, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31260, 3, 3596,
                                                                       3599, 14282, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31270, 3, 3599,
                                                                       3602, 14288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31280, 3, 3602,
                                                                       3605, 14294, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31290, 0, 3,
                                                                       31130, 14204, 31140, 3611,
                                                                       3620, 14300, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31320, 0, 3,
                                                                       31140, 14210, 31150, 3620,
                                                                       3629, 14318, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31350, 0, 3,
                                                                       31150, 14216, 31160, 3629,
                                                                       3638, 14336, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31380, 0, 3,
                                                                       31160, 14222, 31170, 3638,
                                                                       3647, 14354, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31410, 0, 3,
                                                                       31170, 14228, 31180, 3647,
                                                                       3656, 14372, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31440, 0, 3,
                                                                       31180, 14234, 31190, 3656,
                                                                       3665, 14390, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31470, 0, 3,
                                                                       31190, 14240, 31200, 3665,
                                                                       3674, 14408, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31500, 0, 3,
                                                                       31200, 14246, 31210, 3674,
                                                                       3683, 14426, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31530, 0, 3,
                                                                       31210, 14252, 31220, 3683,
                                                                       3692, 14444, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31560, 0, 3,
                                                                       31220, 14258, 31230, 3692,
                                                                       3701, 14462, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31590, 0, 3,
                                                                       31230, 14264, 31240, 3701,
                                                                       3710, 14480, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31620, 0, 3,
                                                                       31240, 14270, 31250, 3710,
                                                                       3719, 14498, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31650, 0, 3,
                                                                       31250, 14276, 31260, 3719,
                                                                       3728, 14516, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31680, 0, 3,
                                                                       31260, 14282, 31270, 3728,
                                                                       3737, 14534, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31710, 0, 3,
                                                                       31270, 14288, 31280, 3737,
                                                                       3746, 14552, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31740, 0, 3,
                                                                       31290, 14300, 31320, 3764,
                                                                       3782, 14570, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31800, 0, 3,
                                                                       31320, 14318, 31350, 3782,
                                                                       3800, 14606, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31860, 0, 3,
                                                                       31350, 14336, 31380, 3800,
                                                                       3818, 14642, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31920, 0, 3,
                                                                       31380, 14354, 31410, 3818,
                                                                       3836, 14678, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31980, 0, 3,
                                                                       31410, 14372, 31440, 3836,
                                                                       3854, 14714, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32040, 0, 3,
                                                                       31440, 14390, 31470, 3854,
                                                                       3872, 14750, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32100, 0, 3,
                                                                       31470, 14408, 31500, 3872,
                                                                       3890, 14786, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32160, 0, 3,
                                                                       31500, 14426, 31530, 3890,
                                                                       3908, 14822, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32220, 0, 3,
                                                                       31530, 14444, 31560, 3908,
                                                                       3926, 14858, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32280, 0, 3,
                                                                       31560, 14462, 31590, 3926,
                                                                       3944, 14894, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32340, 0, 3,
                                                                       31590, 14480, 31620, 3944,
                                                                       3962, 14930, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32400, 0, 3,
                                                                       31620, 14498, 31650, 3962,
                                                                       3980, 14966, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32460, 0, 3,
                                                                       31650, 14516, 31680, 3980,
                                                                       3998, 15002, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32520, 0, 3,
                                                                       31680, 14534, 31710, 3998,
                                                                       4016, 15038, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32580, 0, 3,
                                                                       31740, 14570, 31800, 4052,
                                                                       4082, 15074, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32680, 0, 3,
                                                                       31800, 14606, 31860, 4082,
                                                                       4112, 15134, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32780, 0, 3,
                                                                       31860, 14642, 31920, 4112,
                                                                       4142, 15194, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32880, 0, 3,
                                                                       31920, 14678, 31980, 4142,
                                                                       4172, 15254, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32980, 0, 3,
                                                                       31980, 14714, 32040, 4172,
                                                                       4202, 15314, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33080, 0, 3,
                                                                       32040, 14750, 32100, 4202,
                                                                       4232, 15374, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33180, 0, 3,
                                                                       32100, 14786, 32160, 4232,
                                                                       4262, 15434, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33280, 0, 3,
                                                                       32160, 14822, 32220, 4262,
                                                                       4292, 15494, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33380, 0, 3,
                                                                       32220, 14858, 32280, 4292,
                                                                       4322, 15554, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33480, 0, 3,
                                                                       32280, 14894, 32340, 4322,
                                                                       4352, 15614, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33580, 0, 3,
                                                                       32340, 14930, 32400, 4352,
                                                                       4382, 15674, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33680, 0, 3,
                                                                       32400, 14966, 32460, 4382,
                                                                       4412, 15734, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33780, 0, 3,
                                                                       32460, 15002, 32520, 4412,
                                                                       4442, 15794, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 33880, 0, 3,
                                                                       32580, 15074, 32680, 4502,
                                                                       4547, 15854, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34030, 0, 3,
                                                                       32680, 15134, 32780, 4547,
                                                                       4592, 15944, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34180, 0, 3,
                                                                       32780, 15194, 32880, 4592,
                                                                       4637, 16034, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34330, 0, 3,
                                                                       32880, 15254, 32980, 4637,
                                                                       4682, 16124, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34480, 0, 3,
                                                                       32980, 15314, 33080, 4682,
                                                                       4727, 16214, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34630, 0, 3,
                                                                       33080, 15374, 33180, 4727,
                                                                       4772, 16304, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34780, 0, 3,
                                                                       33180, 15434, 33280, 4772,
                                                                       4817, 16394, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34930, 0, 3,
                                                                       33280, 15494, 33380, 4817,
                                                                       4862, 16484, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35080, 0, 3,
                                                                       33380, 15554, 33480, 4862,
                                                                       4907, 16574, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35230, 0, 3,
                                                                       33480, 15614, 33580, 4907,
                                                                       4952, 16664, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35380, 0, 3,
                                                                       33580, 15674, 33680, 4952,
                                                                       4997, 16754, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35530, 0, 3,
                                                                       33680, 15734, 33780, 4997,
                                                                       5042, 16844, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35680, 0, 3,
                                                                       33880, 15854, 34030, 5132,
                                                                       5195, 16934, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35890, 0, 3,
                                                                       34030, 15944, 34180, 5195,
                                                                       5258, 17060, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36100, 0, 3,
                                                                       34180, 16034, 34330, 5258,
                                                                       5321, 17186, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36310, 0, 3,
                                                                       34330, 16124, 34480, 5321,
                                                                       5384, 17312, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36520, 0, 3,
                                                                       34480, 16214, 34630, 5384,
                                                                       5447, 17438, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36730, 0, 3,
                                                                       34630, 16304, 34780, 5447,
                                                                       5510, 17564, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36940, 0, 3,
                                                                       34780, 16394, 34930, 5510,
                                                                       5573, 17690, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37150, 0, 3,
                                                                       34930, 16484, 35080, 5573,
                                                                       5636, 17816, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37360, 0, 3,
                                                                       35080, 16574, 35230, 5636,
                                                                       5699, 17942, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37570, 0, 3,
                                                                       35230, 16664, 35380, 5699,
                                                                       5762, 18068, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37780, 0, 3,
                                                                       35380, 16754, 35530, 5762,
                                                                       5825, 18194, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37990, 0, 3,
                                                                       35680, 16934, 35890, 5951,
                                                                       6035, 18320, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38270, 0, 3,
                                                                       35890, 17060, 36100, 6035,
                                                                       6119, 18488, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38550, 0, 3,
                                                                       36100, 17186, 36310, 6119,
                                                                       6203, 18656, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38830, 0, 3,
                                                                       36310, 17312, 36520, 6203,
                                                                       6287, 18824, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39110, 0, 3,
                                                                       36520, 17438, 36730, 6287,
                                                                       6371, 18992, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39390, 0, 3,
                                                                       36730, 17564, 36940, 6371,
                                                                       6455, 19160, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39670, 0, 3,
                                                                       36940, 17690, 37150, 6455,
                                                                       6539, 19328, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39950, 0, 3,
                                                                       37150, 17816, 37360, 6539,
                                                                       6623, 19496, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 40230, 0, 3,
                                                                       37360, 17942, 37570, 6623,
                                                                       6707, 19664, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 40510, 0, 3,
                                                                       37570, 18068, 37780, 6707,
                                                                       6791, 19832, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40790, 0, 3,
                                                                       37990, 18320, 38270, 6959,
                                                                       7067, 20000, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41150, 0, 3,
                                                                       38270, 18488, 38550, 7067,
                                                                       7175, 20216, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41510, 0, 3,
                                                                       38550, 18656, 38830, 7175,
                                                                       7283, 20432, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41870, 0, 3,
                                                                       38830, 18824, 39110, 7283,
                                                                       7391, 20648, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42230, 0, 3,
                                                                       39110, 18992, 39390, 7391,
                                                                       7499, 20864, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42590, 0, 3,
                                                                       39390, 19160, 39670, 7499,
                                                                       7607, 21080, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42950, 0, 3,
                                                                       39670, 19328, 39950, 7607,
                                                                       7715, 21296, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 43310, 0, 3,
                                                                       39950, 19496, 40230, 7715,
                                                                       7823, 21512, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 43670, 0, 3,
                                                                       40230, 19664, 40510, 7823,
                                                                       7931, 21728, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44030, 0, 3,
                                                                       40790, 20000, 41150, 8147,
                                                                       8282, 21944, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44480, 0, 3,
                                                                       41150, 20216, 41510, 8282,
                                                                       8417, 22214, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44930, 0, 3,
                                                                       41510, 20432, 41870, 8417,
                                                                       8552, 22484, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 45380, 0, 3,
                                                                       41870, 20648, 42230, 8552,
                                                                       8687, 22754, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 45830, 0, 3,
                                                                       42230, 20864, 42590, 8687,
                                                                       8822, 23024, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 46280, 0, 3,
                                                                       42590, 21080, 42950, 8822,
                                                                       8957, 23294, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 46730, 0, 3,
                                                                       42950, 21296, 43310, 8957,
                                                                       9092, 23564, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 47180, 0, 3,
                                                                       43310, 21512, 43670, 9092,
                                                                       9227, 23834, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 47630, 0, 3,
                                                                       44030, 21944, 44480, 9497,
                                                                       9662, 24104, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 48180, 0, 3,
                                                                       44480, 22214, 44930, 9662,
                                                                       9827, 24434, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 48730, 0, 3,
                                                                       44930, 22484, 45380, 9827,
                                                                       9992, 24764, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 49280, 0, 3,
                                                                       45380, 22754, 45830, 9992,
                                                                       10157, 25094, ncols,
                                                                       gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 49830, 0, 3,
                                                                       45830, 23024, 46280,
                                                                       10157, 10322, 25424,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 50380, 0, 3,
                                                                       46280, 23294, 46730,
                                                                       10322, 10487, 25754,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 50930, 0, 3,
                                                                       46730, 23564, 47180,
                                                                       10487, 10652, 26084,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 51480, 0, 3,
                                                                       47630, 24104, 48180,
                                                                       10982, 11180, 26414,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 52140, 0, 3,
                                                                       48180, 24434, 48730,
                                                                       11180, 11378, 26810,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 52800, 0, 3,
                                                                       48730, 24764, 49280,
                                                                       11378, 11576, 27206,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 53460, 0, 3,
                                                                       49280, 25094, 49830,
                                                                       11576, 11774, 27602,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 54120, 0, 3,
                                                                       49830, 25424, 50380,
                                                                       11774, 11972, 27998,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 54780, 0, 3,
                                                                       50380, 25754, 50930,
                                                                       11972, 12170, 28394,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 55440, 0, 3,
                                                                       51480, 26414, 52140,
                                                                       12566, 12800, 28790,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 56220, 0, 3,
                                                                       52140, 26810, 52800,
                                                                       12800, 13034, 29258,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 57000, 0, 3,
                                                                       52800, 27206, 53460,
                                                                       13034, 13268, 29726,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 57780, 0, 3,
                                                                       53460, 27602, 54120,
                                                                       13268, 13502, 30194,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 58560, 0, 3,
                                                                       54120, 27998, 54780,
                                                                       13502, 13736, 30662,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59340, 3, 14204,
                                                                       14210, 31150, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59355, 3, 14210,
                                                                       14216, 31160, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59370, 3, 14216,
                                                                       14222, 31170, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59385, 3, 14222,
                                                                       14228, 31180, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59400, 3, 14228,
                                                                       14234, 31190, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59415, 3, 14234,
                                                                       14240, 31200, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59430, 3, 14240,
                                                                       14246, 31210, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59445, 3, 14246,
                                                                       14252, 31220, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59460, 3, 14252,
                                                                       14258, 31230, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59475, 3, 14258,
                                                                       14264, 31240, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59490, 3, 14264,
                                                                       14270, 31250, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59505, 3, 14270,
                                                                       14276, 31260, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59520, 3, 14276,
                                                                       14282, 31270, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59535, 3, 14282,
                                                                       14288, 31280, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59550, 0, 3,
                                                                       59340, 31150, 59355,
                                                                       14300, 14318, 31350,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59595, 0, 3,
                                                                       59355, 31160, 59370,
                                                                       14318, 14336, 31380,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59640, 0, 3,
                                                                       59370, 31170, 59385,
                                                                       14336, 14354, 31410,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59685, 0, 3,
                                                                       59385, 31180, 59400,
                                                                       14354, 14372, 31440,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59730, 0, 3,
                                                                       59400, 31190, 59415,
                                                                       14372, 14390, 31470,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59775, 0, 3,
                                                                       59415, 31200, 59430,
                                                                       14390, 14408, 31500,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59820, 0, 3,
                                                                       59430, 31210, 59445,
                                                                       14408, 14426, 31530,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59865, 0, 3,
                                                                       59445, 31220, 59460,
                                                                       14426, 14444, 31560,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59910, 0, 3,
                                                                       59460, 31230, 59475,
                                                                       14444, 14462, 31590,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 59955, 0, 3,
                                                                       59475, 31240, 59490,
                                                                       14462, 14480, 31620,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 60000, 0, 3,
                                                                       59490, 31250, 59505,
                                                                       14480, 14498, 31650,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 60045, 0, 3,
                                                                       59505, 31260, 59520,
                                                                       14498, 14516, 31680,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 60090, 0, 3,
                                                                       59520, 31270, 59535,
                                                                       14516, 14534, 31710,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60135, 0, 3,
                                                                       59550, 31350, 59595,
                                                                       14570, 14606, 31860,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60225, 0, 3,
                                                                       59595, 31380, 59640,
                                                                       14606, 14642, 31920,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60315, 0, 3,
                                                                       59640, 31410, 59685,
                                                                       14642, 14678, 31980,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60405, 0, 3,
                                                                       59685, 31440, 59730,
                                                                       14678, 14714, 32040,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60495, 0, 3,
                                                                       59730, 31470, 59775,
                                                                       14714, 14750, 32100,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60585, 0, 3,
                                                                       59775, 31500, 59820,
                                                                       14750, 14786, 32160,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60675, 0, 3,
                                                                       59820, 31530, 59865,
                                                                       14786, 14822, 32220,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60765, 0, 3,
                                                                       59865, 31560, 59910,
                                                                       14822, 14858, 32280,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60855, 0, 3,
                                                                       59910, 31590, 59955,
                                                                       14858, 14894, 32340,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 60945, 0, 3,
                                                                       59955, 31620, 60000,
                                                                       14894, 14930, 32400,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 61035, 0, 3,
                                                                       60000, 31650, 60045,
                                                                       14930, 14966, 32460,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 61125, 0, 3,
                                                                       60045, 31680, 60090,
                                                                       14966, 15002, 32520,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61215, 0, 3,
                                                                       60135, 31860, 60225,
                                                                       15074, 15134, 32780,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61365, 0, 3,
                                                                       60225, 31920, 60315,
                                                                       15134, 15194, 32880,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61515, 0, 3,
                                                                       60315, 31980, 60405,
                                                                       15194, 15254, 32980,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61665, 0, 3,
                                                                       60405, 32040, 60495,
                                                                       15254, 15314, 33080,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61815, 0, 3,
                                                                       60495, 32100, 60585,
                                                                       15314, 15374, 33180,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61965, 0, 3,
                                                                       60585, 32160, 60675,
                                                                       15374, 15434, 33280,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 62115, 0, 3,
                                                                       60675, 32220, 60765,
                                                                       15434, 15494, 33380,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 62265, 0, 3,
                                                                       60765, 32280, 60855,
                                                                       15494, 15554, 33480,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 62415, 0, 3,
                                                                       60855, 32340, 60945,
                                                                       15554, 15614, 33580,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 62565, 0, 3,
                                                                       60945, 32400, 61035,
                                                                       15614, 15674, 33680,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 62715, 0, 3,
                                                                       61035, 32460, 61125,
                                                                       15674, 15734, 33780,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62865, 0, 3,
                                                                       61215, 32780, 61365,
                                                                       15854, 15944, 34180,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63090, 0, 3,
                                                                       61365, 32880, 61515,
                                                                       15944, 16034, 34330,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63315, 0, 3,
                                                                       61515, 32980, 61665,
                                                                       16034, 16124, 34480,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63540, 0, 3,
                                                                       61665, 33080, 61815,
                                                                       16124, 16214, 34630,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63765, 0, 3,
                                                                       61815, 33180, 61965,
                                                                       16214, 16304, 34780,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63990, 0, 3,
                                                                       61965, 33280, 62115,
                                                                       16304, 16394, 34930,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 64215, 0, 3,
                                                                       62115, 33380, 62265,
                                                                       16394, 16484, 35080,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 64440, 0, 3,
                                                                       62265, 33480, 62415,
                                                                       16484, 16574, 35230,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 64665, 0, 3,
                                                                       62415, 33580, 62565,
                                                                       16574, 16664, 35380,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 64890, 0, 3,
                                                                       62565, 33680, 62715,
                                                                       16664, 16754, 35530,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65115, 0, 3,
                                                                       62865, 34180, 63090,
                                                                       16934, 17060, 36100,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65430, 0, 3,
                                                                       63090, 34330, 63315,
                                                                       17060, 17186, 36310,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65745, 0, 3,
                                                                       63315, 34480, 63540,
                                                                       17186, 17312, 36520,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 66060, 0, 3,
                                                                       63540, 34630, 63765,
                                                                       17312, 17438, 36730,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 66375, 0, 3,
                                                                       63765, 34780, 63990,
                                                                       17438, 17564, 36940,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 66690, 0, 3,
                                                                       63990, 34930, 64215,
                                                                       17564, 17690, 37150,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 67005, 0, 3,
                                                                       64215, 35080, 64440,
                                                                       17690, 17816, 37360,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 67320, 0, 3,
                                                                       64440, 35230, 64665,
                                                                       17816, 17942, 37570,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 67635, 0, 3,
                                                                       64665, 35380, 64890,
                                                                       17942, 18068, 37780,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 67950, 0, 3,
                                                                       65115, 36100, 65430,
                                                                       18320, 18488, 38550,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 68370, 0, 3,
                                                                       65430, 36310, 65745,
                                                                       18488, 18656, 38830,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 68790, 0, 3,
                                                                       65745, 36520, 66060,
                                                                       18656, 18824, 39110,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 69210, 0, 3,
                                                                       66060, 36730, 66375,
                                                                       18824, 18992, 39390,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 69630, 0, 3,
                                                                       66375, 36940, 66690,
                                                                       18992, 19160, 39670,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 70050, 0, 3,
                                                                       66690, 37150, 67005,
                                                                       19160, 19328, 39950,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 70470, 0, 3,
                                                                       67005, 37360, 67320,
                                                                       19328, 19496, 40230,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 70890, 0, 3,
                                                                       67320, 37570, 67635,
                                                                       19496, 19664, 40510,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 71310, 0, 3,
                                                                       67950, 38550, 68370,
                                                                       20000, 20216, 41510,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 71850, 0, 3,
                                                                       68370, 38830, 68790,
                                                                       20216, 20432, 41870,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 72390, 0, 3,
                                                                       68790, 39110, 69210,
                                                                       20432, 20648, 42230,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 72930, 0, 3,
                                                                       69210, 39390, 69630,
                                                                       20648, 20864, 42590,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 73470, 0, 3,
                                                                       69630, 39670, 70050,
                                                                       20864, 21080, 42950,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 74010, 0, 3,
                                                                       70050, 39950, 70470,
                                                                       21080, 21296, 43310,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 74550, 0, 3,
                                                                       70470, 40230, 70890,
                                                                       21296, 21512, 43670,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 75090, 0, 3,
                                                                       71310, 41510, 71850,
                                                                       21944, 22214, 44930,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 75765, 0, 3,
                                                                       71850, 41870, 72390,
                                                                       22214, 22484, 45380,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 76440, 0, 3,
                                                                       72390, 42230, 72930,
                                                                       22484, 22754, 45830,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 77115, 0, 3,
                                                                       72930, 42590, 73470,
                                                                       22754, 23024, 46280,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 77790, 0, 3,
                                                                       73470, 42950, 74010,
                                                                       23024, 23294, 46730,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 78465, 0, 3,
                                                                       74010, 43310, 74550,
                                                                       23294, 23564, 47180,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 79140, 0, 3,
                                                                       75090, 44930, 75765,
                                                                       24104, 24434, 48730,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 79965, 0, 3,
                                                                       75765, 45380, 76440,
                                                                       24434, 24764, 49280,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 80790, 0, 3,
                                                                       76440, 45830, 77115,
                                                                       24764, 25094, 49830,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 81615, 0, 3,
                                                                       77115, 46280, 77790,
                                                                       25094, 25424, 50380,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 82440, 0, 3,
                                                                       77790, 46730, 78465,
                                                                       25424, 25754, 50930,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 83265, 0, 3,
                                                                       79140, 48730, 79965,
                                                                       26414, 26810, 52800,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 84255, 0, 3,
                                                                       79965, 49280, 80790,
                                                                       26810, 27206, 53460,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 85245, 0, 3,
                                                                       80790, 49830, 81615,
                                                                       27206, 27602, 54120,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 86235, 0, 3,
                                                                       81615, 50380, 82440,
                                                                       27602, 27998, 54780,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 87225, 0, 3,
                                                                       83265, 52800, 84255,
                                                                       28790, 29258, 57000,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 88395, 0, 3,
                                                                       84255, 53460, 85245,
                                                                       29258, 29726, 57780,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 89565, 0, 3,
                                                                       85245, 54120, 86235,
                                                                       29726, 30194, 58560,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90735, 3, 31130,
                                                                       31140, 59340, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90756, 3, 31140,
                                                                       31150, 59355, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90777, 3, 31150,
                                                                       31160, 59370, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90798, 3, 31160,
                                                                       31170, 59385, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90819, 3, 31170,
                                                                       31180, 59400, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90840, 3, 31180,
                                                                       31190, 59415, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90861, 3, 31190,
                                                                       31200, 59430, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90882, 3, 31200,
                                                                       31210, 59445, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90903, 3, 31210,
                                                                       31220, 59460, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90924, 3, 31220,
                                                                       31230, 59475, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90945, 3, 31230,
                                                                       31240, 59490, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90966, 3, 31240,
                                                                       31250, 59505, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90987, 3, 31250,
                                                                       31260, 59520, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 91008, 3, 31260,
                                                                       31270, 59535, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91029, 0, 3,
                                                                       90735, 59340, 90756,
                                                                       31290, 31320, 59550,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91092, 0, 3,
                                                                       90756, 59355, 90777,
                                                                       31320, 31350, 59595,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91155, 0, 3,
                                                                       90777, 59370, 90798,
                                                                       31350, 31380, 59640,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91218, 0, 3,
                                                                       90798, 59385, 90819,
                                                                       31380, 31410, 59685,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91281, 0, 3,
                                                                       90819, 59400, 90840,
                                                                       31410, 31440, 59730,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91344, 0, 3,
                                                                       90840, 59415, 90861,
                                                                       31440, 31470, 59775,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91407, 0, 3,
                                                                       90861, 59430, 90882,
                                                                       31470, 31500, 59820,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91470, 0, 3,
                                                                       90882, 59445, 90903,
                                                                       31500, 31530, 59865,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91533, 0, 3,
                                                                       90903, 59460, 90924,
                                                                       31530, 31560, 59910,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91596, 0, 3,
                                                                       90924, 59475, 90945,
                                                                       31560, 31590, 59955,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91659, 0, 3,
                                                                       90945, 59490, 90966,
                                                                       31590, 31620, 60000,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91722, 0, 3,
                                                                       90966, 59505, 90987,
                                                                       31620, 31650, 60045,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 91785, 0, 3,
                                                                       90987, 59520, 91008,
                                                                       31650, 31680, 60090,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91848, 0, 3,
                                                                       91029, 59550, 91092,
                                                                       31740, 31800, 60135,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91974, 0, 3,
                                                                       91092, 59595, 91155,
                                                                       31800, 31860, 60225,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92100, 0, 3,
                                                                       91155, 59640, 91218,
                                                                       31860, 31920, 60315,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92226, 0, 3,
                                                                       91218, 59685, 91281,
                                                                       31920, 31980, 60405,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92352, 0, 3,
                                                                       91281, 59730, 91344,
                                                                       31980, 32040, 60495,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92478, 0, 3,
                                                                       91344, 59775, 91407,
                                                                       32040, 32100, 60585,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92604, 0, 3,
                                                                       91407, 59820, 91470,
                                                                       32100, 32160, 60675,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92730, 0, 3,
                                                                       91470, 59865, 91533,
                                                                       32160, 32220, 60765,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92856, 0, 3,
                                                                       91533, 59910, 91596,
                                                                       32220, 32280, 60855,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 92982, 0, 3,
                                                                       91596, 59955, 91659,
                                                                       32280, 32340, 60945,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 93108, 0, 3,
                                                                       91659, 60000, 91722,
                                                                       32340, 32400, 61035,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 93234, 0, 3,
                                                                       91722, 60045, 91785,
                                                                       32400, 32460, 61125,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93360, 0, 3,
                                                                       91848, 60135, 91974,
                                                                       32580, 32680, 61215,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93570, 0, 3,
                                                                       91974, 60225, 92100,
                                                                       32680, 32780, 61365,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93780, 0, 3,
                                                                       92100, 60315, 92226,
                                                                       32780, 32880, 61515,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93990, 0, 3,
                                                                       92226, 60405, 92352,
                                                                       32880, 32980, 61665,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94200, 0, 3,
                                                                       92352, 60495, 92478,
                                                                       32980, 33080, 61815,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94410, 0, 3,
                                                                       92478, 60585, 92604,
                                                                       33080, 33180, 61965,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94620, 0, 3,
                                                                       92604, 60675, 92730,
                                                                       33180, 33280, 62115,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94830, 0, 3,
                                                                       92730, 60765, 92856,
                                                                       33280, 33380, 62265,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 95040, 0, 3,
                                                                       92856, 60855, 92982,
                                                                       33380, 33480, 62415,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 95250, 0, 3,
                                                                       92982, 60945, 93108,
                                                                       33480, 33580, 62565,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 95460, 0, 3,
                                                                       93108, 61035, 93234,
                                                                       33580, 33680, 62715,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 95670, 0, 3,
                                                                       93360, 61215, 93570,
                                                                       33880, 34030, 62865,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 95985, 0, 3,
                                                                       93570, 61365, 93780,
                                                                       34030, 34180, 63090,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 96300, 0, 3,
                                                                       93780, 61515, 93990,
                                                                       34180, 34330, 63315,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 96615, 0, 3,
                                                                       93990, 61665, 94200,
                                                                       34330, 34480, 63540,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 96930, 0, 3,
                                                                       94200, 61815, 94410,
                                                                       34480, 34630, 63765,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 97245, 0, 3,
                                                                       94410, 61965, 94620,
                                                                       34630, 34780, 63990,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 97560, 0, 3,
                                                                       94620, 62115, 94830,
                                                                       34780, 34930, 64215,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 97875, 0, 3,
                                                                       94830, 62265, 95040,
                                                                       34930, 35080, 64440,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 98190, 0, 3,
                                                                       95040, 62415, 95250,
                                                                       35080, 35230, 64665,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 98505, 0, 3,
                                                                       95250, 62565, 95460,
                                                                       35230, 35380, 64890,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 98820, 0, 3,
                                                                       95670, 62865, 95985,
                                                                       35680, 35890, 65115,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 99261, 0, 3,
                                                                       95985, 63090, 96300,
                                                                       35890, 36100, 65430,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 99702, 0, 3,
                                                                       96300, 63315, 96615,
                                                                       36100, 36310, 65745,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 100143, 0, 3,
                                                                       96615, 63540, 96930,
                                                                       36310, 36520, 66060,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 100584, 0, 3,
                                                                       96930, 63765, 97245,
                                                                       36520, 36730, 66375,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 101025, 0, 3,
                                                                       97245, 63990, 97560,
                                                                       36730, 36940, 66690,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 101466, 0, 3,
                                                                       97560, 64215, 97875,
                                                                       36940, 37150, 67005,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 101907, 0, 3,
                                                                       97875, 64440, 98190,
                                                                       37150, 37360, 67320,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 102348, 0, 3,
                                                                       98190, 64665, 98505,
                                                                       37360, 37570, 67635,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 102789, 0, 3,
                                                                       98820, 65115, 99261,
                                                                       37990, 38270, 67950,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 103377, 0, 3,
                                                                       99261, 65430, 99702,
                                                                       38270, 38550, 68370,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 103965, 0, 3,
                                                                       99702, 65745, 100143,
                                                                       38550, 38830, 68790,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 104553, 0, 3,
                                                                       100143, 66060, 100584,
                                                                       38830, 39110, 69210,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 105141, 0, 3,
                                                                       100584, 66375, 101025,
                                                                       39110, 39390, 69630,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 105729, 0, 3,
                                                                       101025, 66690, 101466,
                                                                       39390, 39670, 70050,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 106317, 0, 3,
                                                                       101466, 67005, 101907,
                                                                       39670, 39950, 70470,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 106905, 0, 3,
                                                                       101907, 67320, 102348,
                                                                       39950, 40230, 70890,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 107493, 0, 3,
                                                                       102789, 67950, 103377,
                                                                       40790, 41150, 71310,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 108249, 0, 3,
                                                                       103377, 68370, 103965,
                                                                       41150, 41510, 71850,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 109005, 0, 3,
                                                                       103965, 68790, 104553,
                                                                       41510, 41870, 72390,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 109761, 0, 3,
                                                                       104553, 69210, 105141,
                                                                       41870, 42230, 72930,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 110517, 0, 3,
                                                                       105141, 69630, 105729,
                                                                       42230, 42590, 73470,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 111273, 0, 3,
                                                                       105729, 70050, 106317,
                                                                       42590, 42950, 74010,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 112029, 0, 3,
                                                                       106317, 70470, 106905,
                                                                       42950, 43310, 74550,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 112785, 0, 3,
                                                                       107493, 71310, 108249,
                                                                       44030, 44480, 75090,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 113730, 0, 3,
                                                                       108249, 71850, 109005,
                                                                       44480, 44930, 75765,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 114675, 0, 3,
                                                                       109005, 72390, 109761,
                                                                       44930, 45380, 76440,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 115620, 0, 3,
                                                                       109761, 72930, 110517,
                                                                       45380, 45830, 77115,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 116565, 0, 3,
                                                                       110517, 73470, 111273,
                                                                       45830, 46280, 77790,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 117510, 0, 3,
                                                                       111273, 74010, 112029,
                                                                       46280, 46730, 78465,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 118455, 0, 3,
                                                                       112785, 75090, 113730,
                                                                       47630, 48180, 79140,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 119610, 0, 3,
                                                                       113730, 75765, 114675,
                                                                       48180, 48730, 79965,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 120765, 0, 3,
                                                                       114675, 76440, 115620,
                                                                       48730, 49280, 80790,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 121920, 0, 3,
                                                                       115620, 77115, 116565,
                                                                       49280, 49830, 81615,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 123075, 0, 3,
                                                                       116565, 77790, 117510,
                                                                       49830, 50380, 82440,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 124230, 0, 3,
                                                                       118455, 79140, 119610,
                                                                       51480, 52140, 83265,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 125616, 0, 3,
                                                                       119610, 79965, 120765,
                                                                       52140, 52800, 84255,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 127002, 0, 3,
                                                                       120765, 80790, 121920,
                                                                       52800, 53460, 85245,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 128388, 0, 3,
                                                                       121920, 81615, 123075,
                                                                       53460, 54120, 86235,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 129774, 0, 3,
                                                                       124230, 83265, 125616,
                                                                       55440, 56220, 87225,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 131412, 0, 3,
                                                                       125616, 84255, 127002,
                                                                       56220, 57000, 88395,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 133050, 0, 3,
                                                                       127002, 85245, 128388,
                                                                       57000, 57780, 89565,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134688, 3, 59340,
                                                                       59355, 90777, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134716, 3, 59355,
                                                                       59370, 90798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134744, 3, 59370,
                                                                       59385, 90819, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134772, 3, 59385,
                                                                       59400, 90840, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134800, 3, 59400,
                                                                       59415, 90861, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134828, 3, 59415,
                                                                       59430, 90882, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134856, 3, 59430,
                                                                       59445, 90903, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134884, 3, 59445,
                                                                       59460, 90924, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134912, 3, 59460,
                                                                       59475, 90945, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134940, 3, 59475,
                                                                       59490, 90966, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134968, 3, 59490,
                                                                       59505, 90987, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134996, 3, 59505,
                                                                       59520, 91008, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135024, 0, 3,
                                                                       134688, 90777, 134716,
                                                                       59550, 59595, 91155,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135108, 0, 3,
                                                                       134716, 90798, 134744,
                                                                       59595, 59640, 91218,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135192, 0, 3,
                                                                       134744, 90819, 134772,
                                                                       59640, 59685, 91281,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135276, 0, 3,
                                                                       134772, 90840, 134800,
                                                                       59685, 59730, 91344,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135360, 0, 3,
                                                                       134800, 90861, 134828,
                                                                       59730, 59775, 91407,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135444, 0, 3,
                                                                       134828, 90882, 134856,
                                                                       59775, 59820, 91470,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135528, 0, 3,
                                                                       134856, 90903, 134884,
                                                                       59820, 59865, 91533,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135612, 0, 3,
                                                                       134884, 90924, 134912,
                                                                       59865, 59910, 91596,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135696, 0, 3,
                                                                       134912, 90945, 134940,
                                                                       59910, 59955, 91659,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135780, 0, 3,
                                                                       134940, 90966, 134968,
                                                                       59955, 60000, 91722,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 135864, 0, 3,
                                                                       134968, 90987, 134996,
                                                                       60000, 60045, 91785,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 135948, 0, 3,
                                                                       135024, 91155, 135108,
                                                                       60135, 60225, 92100,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136116, 0, 3,
                                                                       135108, 91218, 135192,
                                                                       60225, 60315, 92226,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136284, 0, 3,
                                                                       135192, 91281, 135276,
                                                                       60315, 60405, 92352,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136452, 0, 3,
                                                                       135276, 91344, 135360,
                                                                       60405, 60495, 92478,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136620, 0, 3,
                                                                       135360, 91407, 135444,
                                                                       60495, 60585, 92604,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136788, 0, 3,
                                                                       135444, 91470, 135528,
                                                                       60585, 60675, 92730,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 136956, 0, 3,
                                                                       135528, 91533, 135612,
                                                                       60675, 60765, 92856,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 137124, 0, 3,
                                                                       135612, 91596, 135696,
                                                                       60765, 60855, 92982,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 137292, 0, 3,
                                                                       135696, 91659, 135780,
                                                                       60855, 60945, 93108,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 137460, 0, 3,
                                                                       135780, 91722, 135864,
                                                                       60945, 61035, 93234,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 137628, 0, 3,
                                                                       135948, 92100, 136116,
                                                                       61215, 61365, 93780,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 137908, 0, 3,
                                                                       136116, 92226, 136284,
                                                                       61365, 61515, 93990,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 138188, 0, 3,
                                                                       136284, 92352, 136452,
                                                                       61515, 61665, 94200,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 138468, 0, 3,
                                                                       136452, 92478, 136620,
                                                                       61665, 61815, 94410,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 138748, 0, 3,
                                                                       136620, 92604, 136788,
                                                                       61815, 61965, 94620,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 139028, 0, 3,
                                                                       136788, 92730, 136956,
                                                                       61965, 62115, 94830,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 139308, 0, 3,
                                                                       136956, 92856, 137124,
                                                                       62115, 62265, 95040,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 139588, 0, 3,
                                                                       137124, 92982, 137292,
                                                                       62265, 62415, 95250,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 139868, 0, 3,
                                                                       137292, 93108, 137460,
                                                                       62415, 62565, 95460,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 140148, 0, 3,
                                                                       137628, 93780, 137908,
                                                                       62865, 63090, 96300,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 140568, 0, 3,
                                                                       137908, 93990, 138188,
                                                                       63090, 63315, 96615,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 140988, 0, 3,
                                                                       138188, 94200, 138468,
                                                                       63315, 63540, 96930,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 141408, 0, 3,
                                                                       138468, 94410, 138748,
                                                                       63540, 63765, 97245,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 141828, 0, 3,
                                                                       138748, 94620, 139028,
                                                                       63765, 63990, 97560,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 142248, 0, 3,
                                                                       139028, 94830, 139308,
                                                                       63990, 64215, 97875,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 142668, 0, 3,
                                                                       139308, 95040, 139588,
                                                                       64215, 64440, 98190,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 143088, 0, 3,
                                                                       139588, 95250, 139868,
                                                                       64440, 64665, 98505,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 143508, 0, 3,
                                                                       140148, 96300, 140568,
                                                                       65115, 65430, 99702,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 144096, 0, 3,
                                                                       140568, 96615, 140988,
                                                                       65430, 65745, 100143,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 144684, 0, 3,
                                                                       140988, 96930, 141408,
                                                                       65745, 66060, 100584,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 145272, 0, 3,
                                                                       141408, 97245, 141828,
                                                                       66060, 66375, 101025,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 145860, 0, 3,
                                                                       141828, 97560, 142248,
                                                                       66375, 66690, 101466,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 146448, 0, 3,
                                                                       142248, 97875, 142668,
                                                                       66690, 67005, 101907,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 147036, 0, 3,
                                                                       142668, 98190, 143088,
                                                                       67005, 67320, 102348,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 147624, 0, 3,
                                                                       143508, 99702, 144096,
                                                                       67950, 68370, 103965,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 148408, 0, 3,
                                                                       144096, 100143, 144684,
                                                                       68370, 68790, 104553,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 149192, 0, 3,
                                                                       144684, 100584, 145272,
                                                                       68790, 69210, 105141,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 149976, 0, 3,
                                                                       145272, 101025, 145860,
                                                                       69210, 69630, 105729,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 150760, 0, 3,
                                                                       145860, 101466, 146448,
                                                                       69630, 70050, 106317,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 151544, 0, 3,
                                                                       146448, 101907, 147036,
                                                                       70050, 70470, 106905,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 152328, 0, 3,
                                                                       147624, 103965, 148408,
                                                                       71310, 71850, 109005,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 153336, 0, 3,
                                                                       148408, 104553, 149192,
                                                                       71850, 72390, 109761,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 154344, 0, 3,
                                                                       149192, 105141, 149976,
                                                                       72390, 72930, 110517,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 155352, 0, 3,
                                                                       149976, 105729, 150760,
                                                                       72930, 73470, 111273,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 156360, 0, 3,
                                                                       150760, 106317, 151544,
                                                                       73470, 74010, 112029,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 157368, 0, 3,
                                                                       152328, 109005, 153336,
                                                                       75090, 75765, 114675,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 158628, 0, 3,
                                                                       153336, 109761, 154344,
                                                                       75765, 76440, 115620,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 159888, 0, 3,
                                                                       154344, 110517, 155352,
                                                                       76440, 77115, 116565,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 161148, 0, 3,
                                                                       155352, 111273, 156360,
                                                                       77115, 77790, 117510,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 162408, 0, 3,
                                                                       157368, 114675, 158628,
                                                                       79140, 79965, 120765,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 163948, 0, 3,
                                                                       158628, 115620, 159888,
                                                                       79965, 80790, 121920,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 165488, 0, 3,
                                                                       159888, 116565, 161148,
                                                                       80790, 81615, 123075,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 167028, 0, 3,
                                                                       162408, 120765, 163948,
                                                                       83265, 84255, 127002,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 168876, 0, 3,
                                                                       163948, 121920, 165488,
                                                                       84255, 85245, 128388,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 170724, 0, 3,
                                                                       167028, 127002, 168876,
                                                                       87225, 88395, 133050,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172908, 3, 90735,
                                                                       90756, 134688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172944, 3, 90756,
                                                                       90777, 134716, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172980, 3, 90777,
                                                                       90798, 134744, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173016, 3, 90798,
                                                                       90819, 134772, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173052, 3, 90819,
                                                                       90840, 134800, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173088, 3, 90840,
                                                                       90861, 134828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173124, 3, 90861,
                                                                       90882, 134856, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173160, 3, 90882,
                                                                       90903, 134884, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173196, 3, 90903,
                                                                       90924, 134912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173232, 3, 90924,
                                                                       90945, 134940, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173268, 3, 90945,
                                                                       90966, 134968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173304, 3, 90966,
                                                                       90987, 134996, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173340, 0, 3,
                                                                       172908, 134688, 172944,
                                                                       91029, 91092, 135024,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173448, 0, 3,
                                                                       172944, 134716, 172980,
                                                                       91092, 91155, 135108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173556, 0, 3,
                                                                       172980, 134744, 173016,
                                                                       91155, 91218, 135192,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173664, 0, 3,
                                                                       173016, 134772, 173052,
                                                                       91218, 91281, 135276,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173772, 0, 3,
                                                                       173052, 134800, 173088,
                                                                       91281, 91344, 135360,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173880, 0, 3,
                                                                       173088, 134828, 173124,
                                                                       91344, 91407, 135444,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 173988, 0, 3,
                                                                       173124, 134856, 173160,
                                                                       91407, 91470, 135528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 174096, 0, 3,
                                                                       173160, 134884, 173196,
                                                                       91470, 91533, 135612,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 174204, 0, 3,
                                                                       173196, 134912, 173232,
                                                                       91533, 91596, 135696,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 174312, 0, 3,
                                                                       173232, 134940, 173268,
                                                                       91596, 91659, 135780,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 174420, 0, 3,
                                                                       173268, 134968, 173304,
                                                                       91659, 91722, 135864,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 174528, 0, 3,
                                                                       173340, 135024, 173448,
                                                                       91848, 91974, 135948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 174744, 0, 3,
                                                                       173448, 135108, 173556,
                                                                       91974, 92100, 136116,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 174960, 0, 3,
                                                                       173556, 135192, 173664,
                                                                       92100, 92226, 136284,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 175176, 0, 3,
                                                                       173664, 135276, 173772,
                                                                       92226, 92352, 136452,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 175392, 0, 3,
                                                                       173772, 135360, 173880,
                                                                       92352, 92478, 136620,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 175608, 0, 3,
                                                                       173880, 135444, 173988,
                                                                       92478, 92604, 136788,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 175824, 0, 3,
                                                                       173988, 135528, 174096,
                                                                       92604, 92730, 136956,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 176040, 0, 3,
                                                                       174096, 135612, 174204,
                                                                       92730, 92856, 137124,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 176256, 0, 3,
                                                                       174204, 135696, 174312,
                                                                       92856, 92982, 137292,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 176472, 0, 3,
                                                                       174312, 135780, 174420,
                                                                       92982, 93108, 137460,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 176688, 0, 3,
                                                                       174528, 135948, 174744,
                                                                       93360, 93570, 137628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 177048, 0, 3,
                                                                       174744, 136116, 174960,
                                                                       93570, 93780, 137908,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 177408, 0, 3,
                                                                       174960, 136284, 175176,
                                                                       93780, 93990, 138188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 177768, 0, 3,
                                                                       175176, 136452, 175392,
                                                                       93990, 94200, 138468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 178128, 0, 3,
                                                                       175392, 136620, 175608,
                                                                       94200, 94410, 138748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 178488, 0, 3,
                                                                       175608, 136788, 175824,
                                                                       94410, 94620, 139028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 178848, 0, 3,
                                                                       175824, 136956, 176040,
                                                                       94620, 94830, 139308,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 179208, 0, 3,
                                                                       176040, 137124, 176256,
                                                                       94830, 95040, 139588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 179568, 0, 3,
                                                                       176256, 137292, 176472,
                                                                       95040, 95250, 139868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 179928, 0, 3,
                                                                       176688, 137628, 177048,
                                                                       95670, 95985, 140148,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 180468, 0, 3,
                                                                       177048, 137908, 177408,
                                                                       95985, 96300, 140568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 181008, 0, 3,
                                                                       177408, 138188, 177768,
                                                                       96300, 96615, 140988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 181548, 0, 3,
                                                                       177768, 138468, 178128,
                                                                       96615, 96930, 141408,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 182088, 0, 3,
                                                                       178128, 138748, 178488,
                                                                       96930, 97245, 141828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 182628, 0, 3,
                                                                       178488, 139028, 178848,
                                                                       97245, 97560, 142248,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 183168, 0, 3,
                                                                       178848, 139308, 179208,
                                                                       97560, 97875, 142668,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 183708, 0, 3,
                                                                       179208, 139588, 179568,
                                                                       97875, 98190, 143088,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 184248, 0, 3,
                                                                       179928, 140148, 180468,
                                                                       98820, 99261, 143508,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 185004, 0, 3,
                                                                       180468, 140568, 181008,
                                                                       99261, 99702, 144096,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 185760, 0, 3,
                                                                       181008, 140988, 181548,
                                                                       99702, 100143, 144684,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 186516, 0, 3,
                                                                       181548, 141408, 182088,
                                                                       100143, 100584, 145272,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 187272, 0, 3,
                                                                       182088, 141828, 182628,
                                                                       100584, 101025, 145860,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 188028, 0, 3,
                                                                       182628, 142248, 183168,
                                                                       101025, 101466, 146448,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 188784, 0, 3,
                                                                       183168, 142668, 183708,
                                                                       101466, 101907, 147036,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 189540, 0, 3,
                                                                       184248, 143508, 185004,
                                                                       102789, 103377, 147624,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 190548, 0, 3,
                                                                       185004, 144096, 185760,
                                                                       103377, 103965, 148408,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 191556, 0, 3,
                                                                       185760, 144684, 186516,
                                                                       103965, 104553, 149192,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 192564, 0, 3,
                                                                       186516, 145272, 187272,
                                                                       104553, 105141, 149976,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 193572, 0, 3,
                                                                       187272, 145860, 188028,
                                                                       105141, 105729, 150760,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 194580, 0, 3,
                                                                       188028, 146448, 188784,
                                                                       105729, 106317, 151544,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 195588, 0, 3,
                                                                       189540, 147624, 190548,
                                                                       107493, 108249, 152328,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 196884, 0, 3,
                                                                       190548, 148408, 191556,
                                                                       108249, 109005, 153336,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 198180, 0, 3,
                                                                       191556, 149192, 192564,
                                                                       109005, 109761, 154344,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 199476, 0, 3,
                                                                       192564, 149976, 193572,
                                                                       109761, 110517, 155352,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 200772, 0, 3,
                                                                       193572, 150760, 194580,
                                                                       110517, 111273, 156360,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 202068, 0, 3,
                                                                       195588, 152328, 196884,
                                                                       112785, 113730, 157368,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 203688, 0, 3,
                                                                       196884, 153336, 198180,
                                                                       113730, 114675, 158628,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 205308, 0, 3,
                                                                       198180, 154344, 199476,
                                                                       114675, 115620, 159888,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 206928, 0, 3,
                                                                       199476, 155352, 200772,
                                                                       115620, 116565, 161148,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 208548, 0, 3,
                                                                       202068, 157368, 203688,
                                                                       118455, 119610, 162408,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 210528, 0, 3,
                                                                       203688, 158628, 205308,
                                                                       119610, 120765, 163948,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 212508, 0, 3,
                                                                       205308, 159888, 206928,
                                                                       120765, 121920, 165488,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 214488, 0, 3,
                                                                       208548, 162408, 210528,
                                                                       124230, 125616, 167028,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 216864, 0, 3,
                                                                       210528, 163948, 212508,
                                                                       125616, 127002, 168876,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 219240, 0, 3,
                                                                       214488, 167028, 216864,
                                                                       129774, 131412, 170724,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 222048, 189540, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 223476, 195588, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 225312, 202068, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 227607, 208548, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 230412, 214488, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 233778, 219240, 2808, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 223056, 222048, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 224772, 223476, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 226932, 225312, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 229587, 227607, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 232788, 230412, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 236586, 233778, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 237756, 223056, 224772, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 239016, 224772, 226932, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 240636, 226932, 229587, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 242661, 229587, 232788, 15,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 245136, 232788, 236586, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 248106, 237756, 239016, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 250626, 239016, 240636, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 253866, 240636, 242661, 15,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 257916, 242661, 245136, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 262866, 248106, 250626, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 267066, 250626, 253866, 15,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 272466, 253866, 257916, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 279216, 262866, 267066, 15,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 285516, 267066, 272466, 15,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 293616, 279216, 285516, 15,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 302436, 293616, 28, 15, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 302436, 165, nmax);
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
