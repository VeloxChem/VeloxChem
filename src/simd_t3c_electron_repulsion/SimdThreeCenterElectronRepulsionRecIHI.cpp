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


#include "SimdThreeCenterElectronRepulsionRecIHI.hpp"

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

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ihi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ihi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 220844, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1859 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 220844, 148156, 11614, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 17,
                                                             ncols, fj, 6, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3557, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3560, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3563, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3566, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3569, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3572, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3575, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3578, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3581, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3584, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3587, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3590, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3593, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3596, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3599, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3602, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3605, 3, 10, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3614, 3, 11, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3623, 3, 12, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3632, 3, 13, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3641, 3, 14, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3650, 3, 15, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3659, 3, 16, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3668, 3, 17, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3677, 3, 18, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3686, 3, 19, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3695, 3, 20, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3704, 3, 21, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3713, 3, 22, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3722, 3, 23, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3731, 3, 24, 74,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3740, 3, 32, 89,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3758, 3, 35, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3776, 3, 38, 101,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3794, 3, 41, 107,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3812, 3, 44, 113,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3830, 3, 47, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3848, 3, 50, 125,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3866, 3, 53, 131,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3884, 3, 56, 137,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3902, 3, 59, 143,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3920, 3, 62, 149,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3938, 3, 65, 155,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3956, 3, 68, 161,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3974, 3, 71, 167,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3992, 3, 89, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4022, 3, 95, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4052, 3, 101, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4082, 3, 107, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4112, 3, 113, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4142, 3, 119, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4172, 3, 125, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4202, 3, 131, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4232, 3, 137, 273,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4262, 3, 143, 283,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4292, 3, 149, 293,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4322, 3, 155, 303,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4352, 3, 161, 313,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4382, 3, 193, 353,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4427, 3, 203, 368,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4472, 3, 213, 383,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4517, 3, 223, 398,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4562, 3, 233, 413,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4607, 3, 243, 428,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4652, 3, 253, 443,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4697, 3, 263, 458,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4742, 3, 273, 473,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4787, 3, 283, 488,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4832, 3, 293, 503,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4877, 3, 303, 518,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4922, 3, 353, 575,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4985, 3, 368, 596,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5048, 3, 383, 617,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5111, 3, 398, 638,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5174, 3, 413, 659,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5237, 3, 428, 680,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5300, 3, 443, 701,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5363, 3, 458, 722,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5426, 3, 473, 743,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5489, 3, 488, 764,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5552, 3, 503, 785,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5615, 3, 575, 862,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5699, 3, 596, 890,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5783, 3, 617, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5867, 3, 638, 946,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5951, 3, 659, 974,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6035, 3, 680,
                                                                       1002, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6119, 3, 701,
                                                                       1030, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6203, 3, 722,
                                                                       1058, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6287, 3, 743,
                                                                       1086, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6371, 3, 764,
                                                                       1114, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6455, 3, 862,
                                                                       1214, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6563, 3, 890,
                                                                       1250, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6671, 3, 918,
                                                                       1286, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6779, 3, 946,
                                                                       1322, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6887, 3, 974,
                                                                       1358, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6995, 3, 1002,
                                                                       1394, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7103, 3, 1030,
                                                                       1430, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7211, 3, 1058,
                                                                       1466, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7319, 3, 1086,
                                                                       1502, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7427, 3, 1214,
                                                                       1628, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7562, 3, 1250,
                                                                       1673, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7697, 3, 1286,
                                                                       1718, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7832, 3, 1322,
                                                                       1763, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7967, 3, 1358,
                                                                       1808, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8102, 3, 1394,
                                                                       1853, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8237, 3, 1430,
                                                                       1898, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8372, 3, 1466,
                                                                       1943, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8507, 3, 1628,
                                                                       2098, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8672, 3, 1673,
                                                                       2153, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8837, 3, 1718,
                                                                       2208, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9002, 3, 1763,
                                                                       2263, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9167, 3, 1808,
                                                                       2318, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9332, 3, 1853,
                                                                       2373, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9497, 3, 1898,
                                                                       2428, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9662, 3, 2098,
                                                                       2615, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9860, 3, 2153,
                                                                       2681, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10058, 3, 2208,
                                                                       2747, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10256, 3, 2263,
                                                                       2813, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10454, 3, 2318,
                                                                       2879, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10652, 3, 2373,
                                                                       2945, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10850, 3, 2615,
                                                                       3167, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 11084, 3, 2681,
                                                                       3245, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 11318, 3, 2747,
                                                                       3323, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 11552, 3, 2813,
                                                                       3401, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 11786, 3, 2879,
                                                                       3479, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12020, 3, 8, 9,
                                                                       3557, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12026, 3, 9, 10,
                                                                       3560, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12032, 3, 10, 11,
                                                                       3563, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12038, 3, 11, 12,
                                                                       3566, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12044, 3, 12, 13,
                                                                       3569, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12050, 3, 13, 14,
                                                                       3572, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12056, 3, 14, 15,
                                                                       3575, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12062, 3, 15, 16,
                                                                       3578, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12068, 3, 16, 17,
                                                                       3581, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12074, 3, 17, 18,
                                                                       3584, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12080, 3, 18, 19,
                                                                       3587, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12086, 3, 19, 20,
                                                                       3590, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12092, 3, 20, 21,
                                                                       3593, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12098, 3, 21, 22,
                                                                       3596, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12104, 3, 22, 23,
                                                                       3599, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12110, 3, 23, 24,
                                                                       3602, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12116, 0, 3,
                                                                       12020, 3557, 12026, 3605,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12134, 0, 3,
                                                                       12026, 3560, 12032, 3614,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12152, 0, 3,
                                                                       12032, 3563, 12038, 3623,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12170, 0, 3,
                                                                       12038, 3566, 12044, 3632,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12188, 0, 3,
                                                                       12044, 3569, 12050, 3641,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12206, 0, 3,
                                                                       12050, 3572, 12056, 3650,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12224, 0, 3,
                                                                       12056, 3575, 12062, 3659,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12242, 0, 3,
                                                                       12062, 3578, 12068, 3668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12260, 0, 3,
                                                                       12068, 3581, 12074, 3677,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12278, 0, 3,
                                                                       12074, 3584, 12080, 3686,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12296, 0, 3,
                                                                       12080, 3587, 12086, 3695,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12314, 0, 3,
                                                                       12086, 3590, 12092, 3704,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12332, 0, 3,
                                                                       12092, 3593, 12098, 3713,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12350, 0, 3,
                                                                       12098, 3596, 12104, 3722,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12368, 0, 3,
                                                                       12104, 3599, 12110, 3731,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12386, 0, 3,
                                                                       12116, 3605, 12134, 77,
                                                                       83, 3740, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12422, 0, 3,
                                                                       12134, 3614, 12152, 83,
                                                                       89, 3758, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12458, 0, 3,
                                                                       12152, 3623, 12170, 89,
                                                                       95, 3776, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12494, 0, 3,
                                                                       12170, 3632, 12188, 95,
                                                                       101, 3794, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12530, 0, 3,
                                                                       12188, 3641, 12206, 101,
                                                                       107, 3812, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12566, 0, 3,
                                                                       12206, 3650, 12224, 107,
                                                                       113, 3830, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12602, 0, 3,
                                                                       12224, 3659, 12242, 113,
                                                                       119, 3848, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12638, 0, 3,
                                                                       12242, 3668, 12260, 119,
                                                                       125, 3866, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12674, 0, 3,
                                                                       12260, 3677, 12278, 125,
                                                                       131, 3884, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12710, 0, 3,
                                                                       12278, 3686, 12296, 131,
                                                                       137, 3902, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12746, 0, 3,
                                                                       12296, 3695, 12314, 137,
                                                                       143, 3920, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12782, 0, 3,
                                                                       12314, 3704, 12332, 143,
                                                                       149, 3938, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12818, 0, 3,
                                                                       12332, 3713, 12350, 149,
                                                                       155, 3956, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12854, 0, 3,
                                                                       12350, 3722, 12368, 155,
                                                                       161, 3974, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12890, 0, 3,
                                                                       12386, 3740, 12422, 173,
                                                                       183, 3992, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12950, 0, 3,
                                                                       12422, 3758, 12458, 183,
                                                                       193, 4022, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13010, 0, 3,
                                                                       12458, 3776, 12494, 193,
                                                                       203, 4052, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13070, 0, 3,
                                                                       12494, 3794, 12530, 203,
                                                                       213, 4082, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13130, 0, 3,
                                                                       12530, 3812, 12566, 213,
                                                                       223, 4112, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13190, 0, 3,
                                                                       12566, 3830, 12602, 223,
                                                                       233, 4142, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13250, 0, 3,
                                                                       12602, 3848, 12638, 233,
                                                                       243, 4172, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13310, 0, 3,
                                                                       12638, 3866, 12674, 243,
                                                                       253, 4202, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13370, 0, 3,
                                                                       12674, 3884, 12710, 253,
                                                                       263, 4232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13430, 0, 3,
                                                                       12710, 3902, 12746, 263,
                                                                       273, 4262, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13490, 0, 3,
                                                                       12746, 3920, 12782, 273,
                                                                       283, 4292, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13550, 0, 3,
                                                                       12782, 3938, 12818, 283,
                                                                       293, 4322, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13610, 0, 3,
                                                                       12818, 3956, 12854, 293,
                                                                       303, 4352, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13670, 0, 3,
                                                                       12890, 3992, 12950, 323,
                                                                       338, 4382, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13760, 0, 3,
                                                                       12950, 4022, 13010, 338,
                                                                       353, 4427, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13850, 0, 3,
                                                                       13010, 4052, 13070, 353,
                                                                       368, 4472, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13940, 0, 3,
                                                                       13070, 4082, 13130, 368,
                                                                       383, 4517, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14030, 0, 3,
                                                                       13130, 4112, 13190, 383,
                                                                       398, 4562, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14120, 0, 3,
                                                                       13190, 4142, 13250, 398,
                                                                       413, 4607, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14210, 0, 3,
                                                                       13250, 4172, 13310, 413,
                                                                       428, 4652, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14300, 0, 3,
                                                                       13310, 4202, 13370, 428,
                                                                       443, 4697, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14390, 0, 3,
                                                                       13370, 4232, 13430, 443,
                                                                       458, 4742, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14480, 0, 3,
                                                                       13430, 4262, 13490, 458,
                                                                       473, 4787, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14570, 0, 3,
                                                                       13490, 4292, 13550, 473,
                                                                       488, 4832, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14660, 0, 3,
                                                                       13550, 4322, 13610, 488,
                                                                       503, 4877, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14750, 0, 3,
                                                                       13670, 4382, 13760, 533,
                                                                       554, 4922, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14876, 0, 3,
                                                                       13760, 4427, 13850, 554,
                                                                       575, 4985, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15002, 0, 3,
                                                                       13850, 4472, 13940, 575,
                                                                       596, 5048, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15128, 0, 3,
                                                                       13940, 4517, 14030, 596,
                                                                       617, 5111, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15254, 0, 3,
                                                                       14030, 4562, 14120, 617,
                                                                       638, 5174, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15380, 0, 3,
                                                                       14120, 4607, 14210, 638,
                                                                       659, 5237, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15506, 0, 3,
                                                                       14210, 4652, 14300, 659,
                                                                       680, 5300, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15632, 0, 3,
                                                                       14300, 4697, 14390, 680,
                                                                       701, 5363, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15758, 0, 3,
                                                                       14390, 4742, 14480, 701,
                                                                       722, 5426, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15884, 0, 3,
                                                                       14480, 4787, 14570, 722,
                                                                       743, 5489, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16010, 0, 3,
                                                                       14570, 4832, 14660, 743,
                                                                       764, 5552, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16136, 0, 3,
                                                                       14750, 4922, 14876, 806,
                                                                       834, 5615, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16304, 0, 3,
                                                                       14876, 4985, 15002, 834,
                                                                       862, 5699, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16472, 0, 3,
                                                                       15002, 5048, 15128, 862,
                                                                       890, 5783, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16640, 0, 3,
                                                                       15128, 5111, 15254, 890,
                                                                       918, 5867, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16808, 0, 3,
                                                                       15254, 5174, 15380, 918,
                                                                       946, 5951, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16976, 0, 3,
                                                                       15380, 5237, 15506, 946,
                                                                       974, 6035, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17144, 0, 3,
                                                                       15506, 5300, 15632, 974,
                                                                       1002, 6119, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17312, 0, 3,
                                                                       15632, 5363, 15758, 1002,
                                                                       1030, 6203, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17480, 0, 3,
                                                                       15758, 5426, 15884, 1030,
                                                                       1058, 6287, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17648, 0, 3,
                                                                       15884, 5489, 16010, 1058,
                                                                       1086, 6371, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17816, 0, 3,
                                                                       16136, 5615, 16304, 1142,
                                                                       1178, 6455, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18032, 0, 3,
                                                                       16304, 5699, 16472, 1178,
                                                                       1214, 6563, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18248, 0, 3,
                                                                       16472, 5783, 16640, 1214,
                                                                       1250, 6671, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18464, 0, 3,
                                                                       16640, 5867, 16808, 1250,
                                                                       1286, 6779, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18680, 0, 3,
                                                                       16808, 5951, 16976, 1286,
                                                                       1322, 6887, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18896, 0, 3,
                                                                       16976, 6035, 17144, 1322,
                                                                       1358, 6995, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19112, 0, 3,
                                                                       17144, 6119, 17312, 1358,
                                                                       1394, 7103, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19328, 0, 3,
                                                                       17312, 6203, 17480, 1394,
                                                                       1430, 7211, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19544, 0, 3,
                                                                       17480, 6287, 17648, 1430,
                                                                       1466, 7319, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19760, 0, 3,
                                                                       17816, 6455, 18032, 1538,
                                                                       1583, 7427, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20030, 0, 3,
                                                                       18032, 6563, 18248, 1583,
                                                                       1628, 7562, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20300, 0, 3,
                                                                       18248, 6671, 18464, 1628,
                                                                       1673, 7697, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20570, 0, 3,
                                                                       18464, 6779, 18680, 1673,
                                                                       1718, 7832, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20840, 0, 3,
                                                                       18680, 6887, 18896, 1718,
                                                                       1763, 7967, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21110, 0, 3,
                                                                       18896, 6995, 19112, 1763,
                                                                       1808, 8102, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21380, 0, 3,
                                                                       19112, 7103, 19328, 1808,
                                                                       1853, 8237, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21650, 0, 3,
                                                                       19328, 7211, 19544, 1853,
                                                                       1898, 8372, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 21920, 0, 3,
                                                                       19760, 7427, 20030, 1988,
                                                                       2043, 8507, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22250, 0, 3,
                                                                       20030, 7562, 20300, 2043,
                                                                       2098, 8672, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22580, 0, 3,
                                                                       20300, 7697, 20570, 2098,
                                                                       2153, 8837, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22910, 0, 3,
                                                                       20570, 7832, 20840, 2153,
                                                                       2208, 9002, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23240, 0, 3,
                                                                       20840, 7967, 21110, 2208,
                                                                       2263, 9167, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23570, 0, 3,
                                                                       21110, 8102, 21380, 2263,
                                                                       2318, 9332, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23900, 0, 3,
                                                                       21380, 8237, 21650, 2318,
                                                                       2373, 9497, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 24230, 0, 3,
                                                                       21920, 8507, 22250, 2483,
                                                                       2549, 9662, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 24626, 0, 3,
                                                                       22250, 8672, 22580, 2549,
                                                                       2615, 9860, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 25022, 0, 3,
                                                                       22580, 8837, 22910, 2615,
                                                                       2681, 10058, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 25418, 0, 3,
                                                                       22910, 9002, 23240, 2681,
                                                                       2747, 10256, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 25814, 0, 3,
                                                                       23240, 9167, 23570, 2747,
                                                                       2813, 10454, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26210, 0, 3,
                                                                       23570, 9332, 23900, 2813,
                                                                       2879, 10652, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 26606, 0, 3,
                                                                       24230, 9662, 24626, 3011,
                                                                       3089, 10850, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 27074, 0, 3,
                                                                       24626, 9860, 25022, 3089,
                                                                       3167, 11084, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 27542, 0, 3,
                                                                       25022, 10058, 25418, 3167,
                                                                       3245, 11318, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 28010, 0, 3,
                                                                       25418, 10256, 25814, 3245,
                                                                       3323, 11552, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 28478, 0, 3,
                                                                       25814, 10454, 26210, 3323,
                                                                       3401, 11786, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28946, 3, 3557,
                                                                       3560, 12032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28956, 3, 3560,
                                                                       3563, 12038, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28966, 3, 3563,
                                                                       3566, 12044, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28976, 3, 3566,
                                                                       3569, 12050, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28986, 3, 3569,
                                                                       3572, 12056, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28996, 3, 3572,
                                                                       3575, 12062, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29006, 3, 3575,
                                                                       3578, 12068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29016, 3, 3578,
                                                                       3581, 12074, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29026, 3, 3581,
                                                                       3584, 12080, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29036, 3, 3584,
                                                                       3587, 12086, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29046, 3, 3587,
                                                                       3590, 12092, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29056, 3, 3590,
                                                                       3593, 12098, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29066, 3, 3593,
                                                                       3596, 12104, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29076, 3, 3596,
                                                                       3599, 12110, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29086, 0, 3,
                                                                       28946, 12032, 28956, 3605,
                                                                       3614, 12152, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29116, 0, 3,
                                                                       28956, 12038, 28966, 3614,
                                                                       3623, 12170, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29146, 0, 3,
                                                                       28966, 12044, 28976, 3623,
                                                                       3632, 12188, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29176, 0, 3,
                                                                       28976, 12050, 28986, 3632,
                                                                       3641, 12206, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29206, 0, 3,
                                                                       28986, 12056, 28996, 3641,
                                                                       3650, 12224, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29236, 0, 3,
                                                                       28996, 12062, 29006, 3650,
                                                                       3659, 12242, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29266, 0, 3,
                                                                       29006, 12068, 29016, 3659,
                                                                       3668, 12260, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29296, 0, 3,
                                                                       29016, 12074, 29026, 3668,
                                                                       3677, 12278, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29326, 0, 3,
                                                                       29026, 12080, 29036, 3677,
                                                                       3686, 12296, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29356, 0, 3,
                                                                       29036, 12086, 29046, 3686,
                                                                       3695, 12314, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29386, 0, 3,
                                                                       29046, 12092, 29056, 3695,
                                                                       3704, 12332, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29416, 0, 3,
                                                                       29056, 12098, 29066, 3704,
                                                                       3713, 12350, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 29446, 0, 3,
                                                                       29066, 12104, 29076, 3713,
                                                                       3722, 12368, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29476, 0, 3,
                                                                       29086, 12152, 29116, 3740,
                                                                       3758, 12458, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29536, 0, 3,
                                                                       29116, 12170, 29146, 3758,
                                                                       3776, 12494, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29596, 0, 3,
                                                                       29146, 12188, 29176, 3776,
                                                                       3794, 12530, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29656, 0, 3,
                                                                       29176, 12206, 29206, 3794,
                                                                       3812, 12566, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29716, 0, 3,
                                                                       29206, 12224, 29236, 3812,
                                                                       3830, 12602, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29776, 0, 3,
                                                                       29236, 12242, 29266, 3830,
                                                                       3848, 12638, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29836, 0, 3,
                                                                       29266, 12260, 29296, 3848,
                                                                       3866, 12674, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29896, 0, 3,
                                                                       29296, 12278, 29326, 3866,
                                                                       3884, 12710, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29956, 0, 3,
                                                                       29326, 12296, 29356, 3884,
                                                                       3902, 12746, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 30016, 0, 3,
                                                                       29356, 12314, 29386, 3902,
                                                                       3920, 12782, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 30076, 0, 3,
                                                                       29386, 12332, 29416, 3920,
                                                                       3938, 12818, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 30136, 0, 3,
                                                                       29416, 12350, 29446, 3938,
                                                                       3956, 12854, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30196, 0, 3,
                                                                       29476, 12458, 29536, 3992,
                                                                       4022, 13010, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30296, 0, 3,
                                                                       29536, 12494, 29596, 4022,
                                                                       4052, 13070, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30396, 0, 3,
                                                                       29596, 12530, 29656, 4052,
                                                                       4082, 13130, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30496, 0, 3,
                                                                       29656, 12566, 29716, 4082,
                                                                       4112, 13190, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30596, 0, 3,
                                                                       29716, 12602, 29776, 4112,
                                                                       4142, 13250, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30696, 0, 3,
                                                                       29776, 12638, 29836, 4142,
                                                                       4172, 13310, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30796, 0, 3,
                                                                       29836, 12674, 29896, 4172,
                                                                       4202, 13370, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30896, 0, 3,
                                                                       29896, 12710, 29956, 4202,
                                                                       4232, 13430, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30996, 0, 3,
                                                                       29956, 12746, 30016, 4232,
                                                                       4262, 13490, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 31096, 0, 3,
                                                                       30016, 12782, 30076, 4262,
                                                                       4292, 13550, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 31196, 0, 3,
                                                                       30076, 12818, 30136, 4292,
                                                                       4322, 13610, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31296, 0, 3,
                                                                       30196, 13010, 30296, 4382,
                                                                       4427, 13850, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31446, 0, 3,
                                                                       30296, 13070, 30396, 4427,
                                                                       4472, 13940, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31596, 0, 3,
                                                                       30396, 13130, 30496, 4472,
                                                                       4517, 14030, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31746, 0, 3,
                                                                       30496, 13190, 30596, 4517,
                                                                       4562, 14120, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31896, 0, 3,
                                                                       30596, 13250, 30696, 4562,
                                                                       4607, 14210, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32046, 0, 3,
                                                                       30696, 13310, 30796, 4607,
                                                                       4652, 14300, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32196, 0, 3,
                                                                       30796, 13370, 30896, 4652,
                                                                       4697, 14390, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32346, 0, 3,
                                                                       30896, 13430, 30996, 4697,
                                                                       4742, 14480, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32496, 0, 3,
                                                                       30996, 13490, 31096, 4742,
                                                                       4787, 14570, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32646, 0, 3,
                                                                       31096, 13550, 31196, 4787,
                                                                       4832, 14660, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32796, 0, 3,
                                                                       31296, 13850, 31446, 4922,
                                                                       4985, 15002, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33006, 0, 3,
                                                                       31446, 13940, 31596, 4985,
                                                                       5048, 15128, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33216, 0, 3,
                                                                       31596, 14030, 31746, 5048,
                                                                       5111, 15254, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33426, 0, 3,
                                                                       31746, 14120, 31896, 5111,
                                                                       5174, 15380, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33636, 0, 3,
                                                                       31896, 14210, 32046, 5174,
                                                                       5237, 15506, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33846, 0, 3,
                                                                       32046, 14300, 32196, 5237,
                                                                       5300, 15632, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34056, 0, 3,
                                                                       32196, 14390, 32346, 5300,
                                                                       5363, 15758, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34266, 0, 3,
                                                                       32346, 14480, 32496, 5363,
                                                                       5426, 15884, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34476, 0, 3,
                                                                       32496, 14570, 32646, 5426,
                                                                       5489, 16010, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34686, 0, 3,
                                                                       32796, 15002, 33006, 5615,
                                                                       5699, 16472, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34966, 0, 3,
                                                                       33006, 15128, 33216, 5699,
                                                                       5783, 16640, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35246, 0, 3,
                                                                       33216, 15254, 33426, 5783,
                                                                       5867, 16808, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35526, 0, 3,
                                                                       33426, 15380, 33636, 5867,
                                                                       5951, 16976, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35806, 0, 3,
                                                                       33636, 15506, 33846, 5951,
                                                                       6035, 17144, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36086, 0, 3,
                                                                       33846, 15632, 34056, 6035,
                                                                       6119, 17312, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36366, 0, 3,
                                                                       34056, 15758, 34266, 6119,
                                                                       6203, 17480, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36646, 0, 3,
                                                                       34266, 15884, 34476, 6203,
                                                                       6287, 17648, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36926, 0, 3,
                                                                       34686, 16472, 34966, 6455,
                                                                       6563, 18248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 37286, 0, 3,
                                                                       34966, 16640, 35246, 6563,
                                                                       6671, 18464, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 37646, 0, 3,
                                                                       35246, 16808, 35526, 6671,
                                                                       6779, 18680, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 38006, 0, 3,
                                                                       35526, 16976, 35806, 6779,
                                                                       6887, 18896, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 38366, 0, 3,
                                                                       35806, 17144, 36086, 6887,
                                                                       6995, 19112, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 38726, 0, 3,
                                                                       36086, 17312, 36366, 6995,
                                                                       7103, 19328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 39086, 0, 3,
                                                                       36366, 17480, 36646, 7103,
                                                                       7211, 19544, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 39446, 0, 3,
                                                                       36926, 18248, 37286, 7427,
                                                                       7562, 20300, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 39896, 0, 3,
                                                                       37286, 18464, 37646, 7562,
                                                                       7697, 20570, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 40346, 0, 3,
                                                                       37646, 18680, 38006, 7697,
                                                                       7832, 20840, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 40796, 0, 3,
                                                                       38006, 18896, 38366, 7832,
                                                                       7967, 21110, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 41246, 0, 3,
                                                                       38366, 19112, 38726, 7967,
                                                                       8102, 21380, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 41696, 0, 3,
                                                                       38726, 19328, 39086, 8102,
                                                                       8237, 21650, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 42146, 0, 3,
                                                                       39446, 20300, 39896, 8507,
                                                                       8672, 22580, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 42696, 0, 3,
                                                                       39896, 20570, 40346, 8672,
                                                                       8837, 22910, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 43246, 0, 3,
                                                                       40346, 20840, 40796, 8837,
                                                                       9002, 23240, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 43796, 0, 3,
                                                                       40796, 21110, 41246, 9002,
                                                                       9167, 23570, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 44346, 0, 3,
                                                                       41246, 21380, 41696, 9167,
                                                                       9332, 23900, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 44896, 0, 3,
                                                                       42146, 22580, 42696, 9662,
                                                                       9860, 25022, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 45556, 0, 3,
                                                                       42696, 22910, 43246, 9860,
                                                                       10058, 25418, ncols,
                                                                       gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 46216, 0, 3,
                                                                       43246, 23240, 43796,
                                                                       10058, 10256, 25814,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 46876, 0, 3,
                                                                       43796, 23570, 44346,
                                                                       10256, 10454, 26210,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 47536, 0, 3,
                                                                       44896, 25022, 45556,
                                                                       10850, 11084, 27542,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 48316, 0, 3,
                                                                       45556, 25418, 46216,
                                                                       11084, 11318, 28010,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 49096, 0, 3,
                                                                       46216, 25814, 46876,
                                                                       11318, 11552, 28478,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49876, 3, 12020,
                                                                       12026, 28946, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49891, 3, 12026,
                                                                       12032, 28956, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49906, 3, 12032,
                                                                       12038, 28966, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49921, 3, 12038,
                                                                       12044, 28976, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49936, 3, 12044,
                                                                       12050, 28986, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49951, 3, 12050,
                                                                       12056, 28996, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49966, 3, 12056,
                                                                       12062, 29006, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49981, 3, 12062,
                                                                       12068, 29016, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49996, 3, 12068,
                                                                       12074, 29026, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 50011, 3, 12074,
                                                                       12080, 29036, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 50026, 3, 12080,
                                                                       12086, 29046, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 50041, 3, 12086,
                                                                       12092, 29056, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 50056, 3, 12092,
                                                                       12098, 29066, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 50071, 3, 12098,
                                                                       12104, 29076, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50086, 0, 3,
                                                                       49876, 28946, 49891,
                                                                       12116, 12134, 29086,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50131, 0, 3,
                                                                       49891, 28956, 49906,
                                                                       12134, 12152, 29116,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50176, 0, 3,
                                                                       49906, 28966, 49921,
                                                                       12152, 12170, 29146,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50221, 0, 3,
                                                                       49921, 28976, 49936,
                                                                       12170, 12188, 29176,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50266, 0, 3,
                                                                       49936, 28986, 49951,
                                                                       12188, 12206, 29206,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50311, 0, 3,
                                                                       49951, 28996, 49966,
                                                                       12206, 12224, 29236,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50356, 0, 3,
                                                                       49966, 29006, 49981,
                                                                       12224, 12242, 29266,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50401, 0, 3,
                                                                       49981, 29016, 49996,
                                                                       12242, 12260, 29296,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50446, 0, 3,
                                                                       49996, 29026, 50011,
                                                                       12260, 12278, 29326,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50491, 0, 3,
                                                                       50011, 29036, 50026,
                                                                       12278, 12296, 29356,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50536, 0, 3,
                                                                       50026, 29046, 50041,
                                                                       12296, 12314, 29386,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50581, 0, 3,
                                                                       50041, 29056, 50056,
                                                                       12314, 12332, 29416,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50626, 0, 3,
                                                                       50056, 29066, 50071,
                                                                       12332, 12350, 29446,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50671, 0, 3,
                                                                       50086, 29086, 50131,
                                                                       12386, 12422, 29476,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50761, 0, 3,
                                                                       50131, 29116, 50176,
                                                                       12422, 12458, 29536,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50851, 0, 3,
                                                                       50176, 29146, 50221,
                                                                       12458, 12494, 29596,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50941, 0, 3,
                                                                       50221, 29176, 50266,
                                                                       12494, 12530, 29656,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51031, 0, 3,
                                                                       50266, 29206, 50311,
                                                                       12530, 12566, 29716,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51121, 0, 3,
                                                                       50311, 29236, 50356,
                                                                       12566, 12602, 29776,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51211, 0, 3,
                                                                       50356, 29266, 50401,
                                                                       12602, 12638, 29836,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51301, 0, 3,
                                                                       50401, 29296, 50446,
                                                                       12638, 12674, 29896,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51391, 0, 3,
                                                                       50446, 29326, 50491,
                                                                       12674, 12710, 29956,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51481, 0, 3,
                                                                       50491, 29356, 50536,
                                                                       12710, 12746, 30016,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51571, 0, 3,
                                                                       50536, 29386, 50581,
                                                                       12746, 12782, 30076,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51661, 0, 3,
                                                                       50581, 29416, 50626,
                                                                       12782, 12818, 30136,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51751, 0, 3,
                                                                       50671, 29476, 50761,
                                                                       12890, 12950, 30196,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51901, 0, 3,
                                                                       50761, 29536, 50851,
                                                                       12950, 13010, 30296,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52051, 0, 3,
                                                                       50851, 29596, 50941,
                                                                       13010, 13070, 30396,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52201, 0, 3,
                                                                       50941, 29656, 51031,
                                                                       13070, 13130, 30496,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52351, 0, 3,
                                                                       51031, 29716, 51121,
                                                                       13130, 13190, 30596,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52501, 0, 3,
                                                                       51121, 29776, 51211,
                                                                       13190, 13250, 30696,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52651, 0, 3,
                                                                       51211, 29836, 51301,
                                                                       13250, 13310, 30796,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52801, 0, 3,
                                                                       51301, 29896, 51391,
                                                                       13310, 13370, 30896,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52951, 0, 3,
                                                                       51391, 29956, 51481,
                                                                       13370, 13430, 30996,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53101, 0, 3,
                                                                       51481, 30016, 51571,
                                                                       13430, 13490, 31096,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53251, 0, 3,
                                                                       51571, 30076, 51661,
                                                                       13490, 13550, 31196,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53401, 0, 3,
                                                                       51751, 30196, 51901,
                                                                       13670, 13760, 31296,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53626, 0, 3,
                                                                       51901, 30296, 52051,
                                                                       13760, 13850, 31446,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53851, 0, 3,
                                                                       52051, 30396, 52201,
                                                                       13850, 13940, 31596,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54076, 0, 3,
                                                                       52201, 30496, 52351,
                                                                       13940, 14030, 31746,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54301, 0, 3,
                                                                       52351, 30596, 52501,
                                                                       14030, 14120, 31896,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54526, 0, 3,
                                                                       52501, 30696, 52651,
                                                                       14120, 14210, 32046,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54751, 0, 3,
                                                                       52651, 30796, 52801,
                                                                       14210, 14300, 32196,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54976, 0, 3,
                                                                       52801, 30896, 52951,
                                                                       14300, 14390, 32346,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55201, 0, 3,
                                                                       52951, 30996, 53101,
                                                                       14390, 14480, 32496,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55426, 0, 3,
                                                                       53101, 31096, 53251,
                                                                       14480, 14570, 32646,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 55651, 0, 3,
                                                                       53401, 31296, 53626,
                                                                       14750, 14876, 32796,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 55966, 0, 3,
                                                                       53626, 31446, 53851,
                                                                       14876, 15002, 33006,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 56281, 0, 3,
                                                                       53851, 31596, 54076,
                                                                       15002, 15128, 33216,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 56596, 0, 3,
                                                                       54076, 31746, 54301,
                                                                       15128, 15254, 33426,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 56911, 0, 3,
                                                                       54301, 31896, 54526,
                                                                       15254, 15380, 33636,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57226, 0, 3,
                                                                       54526, 32046, 54751,
                                                                       15380, 15506, 33846,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57541, 0, 3,
                                                                       54751, 32196, 54976,
                                                                       15506, 15632, 34056,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57856, 0, 3,
                                                                       54976, 32346, 55201,
                                                                       15632, 15758, 34266,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58171, 0, 3,
                                                                       55201, 32496, 55426,
                                                                       15758, 15884, 34476,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 58486, 0, 3,
                                                                       55651, 32796, 55966,
                                                                       16136, 16304, 34686,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 58906, 0, 3,
                                                                       55966, 33006, 56281,
                                                                       16304, 16472, 34966,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 59326, 0, 3,
                                                                       56281, 33216, 56596,
                                                                       16472, 16640, 35246,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 59746, 0, 3,
                                                                       56596, 33426, 56911,
                                                                       16640, 16808, 35526,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 60166, 0, 3,
                                                                       56911, 33636, 57226,
                                                                       16808, 16976, 35806,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 60586, 0, 3,
                                                                       57226, 33846, 57541,
                                                                       16976, 17144, 36086,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61006, 0, 3,
                                                                       57541, 34056, 57856,
                                                                       17144, 17312, 36366,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61426, 0, 3,
                                                                       57856, 34266, 58171,
                                                                       17312, 17480, 36646,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 61846, 0, 3,
                                                                       58486, 34686, 58906,
                                                                       17816, 18032, 36926,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 62386, 0, 3,
                                                                       58906, 34966, 59326,
                                                                       18032, 18248, 37286,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 62926, 0, 3,
                                                                       59326, 35246, 59746,
                                                                       18248, 18464, 37646,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 63466, 0, 3,
                                                                       59746, 35526, 60166,
                                                                       18464, 18680, 38006,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 64006, 0, 3,
                                                                       60166, 35806, 60586,
                                                                       18680, 18896, 38366,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 64546, 0, 3,
                                                                       60586, 36086, 61006,
                                                                       18896, 19112, 38726,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 65086, 0, 3,
                                                                       61006, 36366, 61426,
                                                                       19112, 19328, 39086,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 65626, 0, 3,
                                                                       61846, 36926, 62386,
                                                                       19760, 20030, 39446,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 66301, 0, 3,
                                                                       62386, 37286, 62926,
                                                                       20030, 20300, 39896,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 66976, 0, 3,
                                                                       62926, 37646, 63466,
                                                                       20300, 20570, 40346,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 67651, 0, 3,
                                                                       63466, 38006, 64006,
                                                                       20570, 20840, 40796,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 68326, 0, 3,
                                                                       64006, 38366, 64546,
                                                                       20840, 21110, 41246,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 69001, 0, 3,
                                                                       64546, 38726, 65086,
                                                                       21110, 21380, 41696,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 69676, 0, 3,
                                                                       65626, 39446, 66301,
                                                                       21920, 22250, 42146,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 70501, 0, 3,
                                                                       66301, 39896, 66976,
                                                                       22250, 22580, 42696,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 71326, 0, 3,
                                                                       66976, 40346, 67651,
                                                                       22580, 22910, 43246,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 72151, 0, 3,
                                                                       67651, 40796, 68326,
                                                                       22910, 23240, 43796,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 72976, 0, 3,
                                                                       68326, 41246, 69001,
                                                                       23240, 23570, 44346,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 73801, 0, 3,
                                                                       69676, 42146, 70501,
                                                                       24230, 24626, 44896,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 74791, 0, 3,
                                                                       70501, 42696, 71326,
                                                                       24626, 25022, 45556,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 75781, 0, 3,
                                                                       71326, 43246, 72151,
                                                                       25022, 25418, 46216,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 76771, 0, 3,
                                                                       72151, 43796, 72976,
                                                                       25418, 25814, 46876,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 77761, 0, 3,
                                                                       73801, 44896, 74791,
                                                                       26606, 27074, 47536,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 78931, 0, 3,
                                                                       74791, 45556, 75781,
                                                                       27074, 27542, 48316,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 80101, 0, 3,
                                                                       75781, 46216, 76771,
                                                                       27542, 28010, 49096,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81271, 3, 28946,
                                                                       28956, 49906, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81292, 3, 28956,
                                                                       28966, 49921, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81313, 3, 28966,
                                                                       28976, 49936, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81334, 3, 28976,
                                                                       28986, 49951, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81355, 3, 28986,
                                                                       28996, 49966, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81376, 3, 28996,
                                                                       29006, 49981, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81397, 3, 29006,
                                                                       29016, 49996, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81418, 3, 29016,
                                                                       29026, 50011, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81439, 3, 29026,
                                                                       29036, 50026, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81460, 3, 29036,
                                                                       29046, 50041, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81481, 3, 29046,
                                                                       29056, 50056, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81502, 3, 29056,
                                                                       29066, 50071, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81523, 0, 3,
                                                                       81271, 49906, 81292,
                                                                       29086, 29116, 50176,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81586, 0, 3,
                                                                       81292, 49921, 81313,
                                                                       29116, 29146, 50221,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81649, 0, 3,
                                                                       81313, 49936, 81334,
                                                                       29146, 29176, 50266,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81712, 0, 3,
                                                                       81334, 49951, 81355,
                                                                       29176, 29206, 50311,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81775, 0, 3,
                                                                       81355, 49966, 81376,
                                                                       29206, 29236, 50356,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81838, 0, 3,
                                                                       81376, 49981, 81397,
                                                                       29236, 29266, 50401,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81901, 0, 3,
                                                                       81397, 49996, 81418,
                                                                       29266, 29296, 50446,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81964, 0, 3,
                                                                       81418, 50011, 81439,
                                                                       29296, 29326, 50491,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82027, 0, 3,
                                                                       81439, 50026, 81460,
                                                                       29326, 29356, 50536,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82090, 0, 3,
                                                                       81460, 50041, 81481,
                                                                       29356, 29386, 50581,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82153, 0, 3,
                                                                       81481, 50056, 81502,
                                                                       29386, 29416, 50626,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82216, 0, 3,
                                                                       81523, 50176, 81586,
                                                                       29476, 29536, 50851,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82342, 0, 3,
                                                                       81586, 50221, 81649,
                                                                       29536, 29596, 50941,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82468, 0, 3,
                                                                       81649, 50266, 81712,
                                                                       29596, 29656, 51031,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82594, 0, 3,
                                                                       81712, 50311, 81775,
                                                                       29656, 29716, 51121,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82720, 0, 3,
                                                                       81775, 50356, 81838,
                                                                       29716, 29776, 51211,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82846, 0, 3,
                                                                       81838, 50401, 81901,
                                                                       29776, 29836, 51301,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82972, 0, 3,
                                                                       81901, 50446, 81964,
                                                                       29836, 29896, 51391,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83098, 0, 3,
                                                                       81964, 50491, 82027,
                                                                       29896, 29956, 51481,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83224, 0, 3,
                                                                       82027, 50536, 82090,
                                                                       29956, 30016, 51571,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83350, 0, 3,
                                                                       82090, 50581, 82153,
                                                                       30016, 30076, 51661,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 83476, 0, 3,
                                                                       82216, 50851, 82342,
                                                                       30196, 30296, 52051,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 83686, 0, 3,
                                                                       82342, 50941, 82468,
                                                                       30296, 30396, 52201,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 83896, 0, 3,
                                                                       82468, 51031, 82594,
                                                                       30396, 30496, 52351,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84106, 0, 3,
                                                                       82594, 51121, 82720,
                                                                       30496, 30596, 52501,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84316, 0, 3,
                                                                       82720, 51211, 82846,
                                                                       30596, 30696, 52651,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84526, 0, 3,
                                                                       82846, 51301, 82972,
                                                                       30696, 30796, 52801,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84736, 0, 3,
                                                                       82972, 51391, 83098,
                                                                       30796, 30896, 52951,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84946, 0, 3,
                                                                       83098, 51481, 83224,
                                                                       30896, 30996, 53101,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 85156, 0, 3,
                                                                       83224, 51571, 83350,
                                                                       30996, 31096, 53251,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 85366, 0, 3,
                                                                       83476, 52051, 83686,
                                                                       31296, 31446, 53851,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 85681, 0, 3,
                                                                       83686, 52201, 83896,
                                                                       31446, 31596, 54076,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 85996, 0, 3,
                                                                       83896, 52351, 84106,
                                                                       31596, 31746, 54301,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 86311, 0, 3,
                                                                       84106, 52501, 84316,
                                                                       31746, 31896, 54526,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 86626, 0, 3,
                                                                       84316, 52651, 84526,
                                                                       31896, 32046, 54751,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 86941, 0, 3,
                                                                       84526, 52801, 84736,
                                                                       32046, 32196, 54976,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 87256, 0, 3,
                                                                       84736, 52951, 84946,
                                                                       32196, 32346, 55201,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 87571, 0, 3,
                                                                       84946, 53101, 85156,
                                                                       32346, 32496, 55426,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 87886, 0, 3,
                                                                       85366, 53851, 85681,
                                                                       32796, 33006, 56281,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 88327, 0, 3,
                                                                       85681, 54076, 85996,
                                                                       33006, 33216, 56596,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 88768, 0, 3,
                                                                       85996, 54301, 86311,
                                                                       33216, 33426, 56911,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 89209, 0, 3,
                                                                       86311, 54526, 86626,
                                                                       33426, 33636, 57226,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 89650, 0, 3,
                                                                       86626, 54751, 86941,
                                                                       33636, 33846, 57541,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 90091, 0, 3,
                                                                       86941, 54976, 87256,
                                                                       33846, 34056, 57856,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 90532, 0, 3,
                                                                       87256, 55201, 87571,
                                                                       34056, 34266, 58171,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 90973, 0, 3,
                                                                       87886, 56281, 88327,
                                                                       34686, 34966, 59326,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 91561, 0, 3,
                                                                       88327, 56596, 88768,
                                                                       34966, 35246, 59746,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 92149, 0, 3,
                                                                       88768, 56911, 89209,
                                                                       35246, 35526, 60166,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 92737, 0, 3,
                                                                       89209, 57226, 89650,
                                                                       35526, 35806, 60586,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 93325, 0, 3,
                                                                       89650, 57541, 90091,
                                                                       35806, 36086, 61006,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 93913, 0, 3,
                                                                       90091, 57856, 90532,
                                                                       36086, 36366, 61426,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 94501, 0, 3,
                                                                       90973, 59326, 91561,
                                                                       36926, 37286, 62926,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 95257, 0, 3,
                                                                       91561, 59746, 92149,
                                                                       37286, 37646, 63466,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 96013, 0, 3,
                                                                       92149, 60166, 92737,
                                                                       37646, 38006, 64006,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 96769, 0, 3,
                                                                       92737, 60586, 93325,
                                                                       38006, 38366, 64546,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 97525, 0, 3,
                                                                       93325, 61006, 93913,
                                                                       38366, 38726, 65086,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 98281, 0, 3,
                                                                       94501, 62926, 95257,
                                                                       39446, 39896, 66976,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 99226, 0, 3,
                                                                       95257, 63466, 96013,
                                                                       39896, 40346, 67651,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 100171, 0, 3,
                                                                       96013, 64006, 96769,
                                                                       40346, 40796, 68326,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 101116, 0, 3,
                                                                       96769, 64546, 97525,
                                                                       40796, 41246, 69001,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 102061, 0, 3,
                                                                       98281, 66976, 99226,
                                                                       42146, 42696, 71326,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 103216, 0, 3,
                                                                       99226, 67651, 100171,
                                                                       42696, 43246, 72151,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 104371, 0, 3,
                                                                       100171, 68326, 101116,
                                                                       43246, 43796, 72976,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 105526, 0, 3,
                                                                       102061, 71326, 103216,
                                                                       44896, 45556, 75781,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 106912, 0, 3,
                                                                       103216, 72151, 104371,
                                                                       45556, 46216, 76771,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 108298, 0, 3,
                                                                       105526, 75781, 106912,
                                                                       47536, 48316, 80101,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 109936, 3, 49876,
                                                                       49891, 81271, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 109964, 3, 49891,
                                                                       49906, 81292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 109992, 3, 49906,
                                                                       49921, 81313, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110020, 3, 49921,
                                                                       49936, 81334, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110048, 3, 49936,
                                                                       49951, 81355, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110076, 3, 49951,
                                                                       49966, 81376, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110104, 3, 49966,
                                                                       49981, 81397, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110132, 3, 49981,
                                                                       49996, 81418, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110160, 3, 49996,
                                                                       50011, 81439, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110188, 3, 50011,
                                                                       50026, 81460, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110216, 3, 50026,
                                                                       50041, 81481, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110244, 3, 50041,
                                                                       50056, 81502, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 110272, 0, 3,
                                                                       109936, 81271, 109964,
                                                                       50086, 50131, 81523,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 110356, 0, 3,
                                                                       109964, 81292, 109992,
                                                                       50131, 50176, 81586,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 110440, 0, 3,
                                                                       109992, 81313, 110020,
                                                                       50176, 50221, 81649,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 110524, 0, 3,
                                                                       110020, 81334, 110048,
                                                                       50221, 50266, 81712,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 110608, 0, 3,
                                                                       110048, 81355, 110076,
                                                                       50266, 50311, 81775,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 110692, 0, 3,
                                                                       110076, 81376, 110104,
                                                                       50311, 50356, 81838,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 110776, 0, 3,
                                                                       110104, 81397, 110132,
                                                                       50356, 50401, 81901,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 110860, 0, 3,
                                                                       110132, 81418, 110160,
                                                                       50401, 50446, 81964,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 110944, 0, 3,
                                                                       110160, 81439, 110188,
                                                                       50446, 50491, 82027,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 111028, 0, 3,
                                                                       110188, 81460, 110216,
                                                                       50491, 50536, 82090,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 111112, 0, 3,
                                                                       110216, 81481, 110244,
                                                                       50536, 50581, 82153,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 111196, 0, 3,
                                                                       110272, 81523, 110356,
                                                                       50671, 50761, 82216,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 111364, 0, 3,
                                                                       110356, 81586, 110440,
                                                                       50761, 50851, 82342,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 111532, 0, 3,
                                                                       110440, 81649, 110524,
                                                                       50851, 50941, 82468,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 111700, 0, 3,
                                                                       110524, 81712, 110608,
                                                                       50941, 51031, 82594,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 111868, 0, 3,
                                                                       110608, 81775, 110692,
                                                                       51031, 51121, 82720,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 112036, 0, 3,
                                                                       110692, 81838, 110776,
                                                                       51121, 51211, 82846,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 112204, 0, 3,
                                                                       110776, 81901, 110860,
                                                                       51211, 51301, 82972,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 112372, 0, 3,
                                                                       110860, 81964, 110944,
                                                                       51301, 51391, 83098,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 112540, 0, 3,
                                                                       110944, 82027, 111028,
                                                                       51391, 51481, 83224,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 112708, 0, 3,
                                                                       111028, 82090, 111112,
                                                                       51481, 51571, 83350,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 112876, 0, 3,
                                                                       111196, 82216, 111364,
                                                                       51751, 51901, 83476,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 113156, 0, 3,
                                                                       111364, 82342, 111532,
                                                                       51901, 52051, 83686,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 113436, 0, 3,
                                                                       111532, 82468, 111700,
                                                                       52051, 52201, 83896,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 113716, 0, 3,
                                                                       111700, 82594, 111868,
                                                                       52201, 52351, 84106,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 113996, 0, 3,
                                                                       111868, 82720, 112036,
                                                                       52351, 52501, 84316,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 114276, 0, 3,
                                                                       112036, 82846, 112204,
                                                                       52501, 52651, 84526,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 114556, 0, 3,
                                                                       112204, 82972, 112372,
                                                                       52651, 52801, 84736,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 114836, 0, 3,
                                                                       112372, 83098, 112540,
                                                                       52801, 52951, 84946,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 115116, 0, 3,
                                                                       112540, 83224, 112708,
                                                                       52951, 53101, 85156,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 115396, 0, 3,
                                                                       112876, 83476, 113156,
                                                                       53401, 53626, 85366,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 115816, 0, 3,
                                                                       113156, 83686, 113436,
                                                                       53626, 53851, 85681,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 116236, 0, 3,
                                                                       113436, 83896, 113716,
                                                                       53851, 54076, 85996,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 116656, 0, 3,
                                                                       113716, 84106, 113996,
                                                                       54076, 54301, 86311,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 117076, 0, 3,
                                                                       113996, 84316, 114276,
                                                                       54301, 54526, 86626,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 117496, 0, 3,
                                                                       114276, 84526, 114556,
                                                                       54526, 54751, 86941,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 117916, 0, 3,
                                                                       114556, 84736, 114836,
                                                                       54751, 54976, 87256,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 118336, 0, 3,
                                                                       114836, 84946, 115116,
                                                                       54976, 55201, 87571,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 118756, 0, 3,
                                                                       115396, 85366, 115816,
                                                                       55651, 55966, 87886,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 119344, 0, 3,
                                                                       115816, 85681, 116236,
                                                                       55966, 56281, 88327,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 119932, 0, 3,
                                                                       116236, 85996, 116656,
                                                                       56281, 56596, 88768,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 120520, 0, 3,
                                                                       116656, 86311, 117076,
                                                                       56596, 56911, 89209,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 121108, 0, 3,
                                                                       117076, 86626, 117496,
                                                                       56911, 57226, 89650,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 121696, 0, 3,
                                                                       117496, 86941, 117916,
                                                                       57226, 57541, 90091,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 122284, 0, 3,
                                                                       117916, 87256, 118336,
                                                                       57541, 57856, 90532,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 122872, 0, 3,
                                                                       118756, 87886, 119344,
                                                                       58486, 58906, 90973,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 123656, 0, 3,
                                                                       119344, 88327, 119932,
                                                                       58906, 59326, 91561,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 124440, 0, 3,
                                                                       119932, 88768, 120520,
                                                                       59326, 59746, 92149,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 125224, 0, 3,
                                                                       120520, 89209, 121108,
                                                                       59746, 60166, 92737,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 126008, 0, 3,
                                                                       121108, 89650, 121696,
                                                                       60166, 60586, 93325,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 126792, 0, 3,
                                                                       121696, 90091, 122284,
                                                                       60586, 61006, 93913,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 127576, 0, 3,
                                                                       122872, 90973, 123656,
                                                                       61846, 62386, 94501,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 128584, 0, 3,
                                                                       123656, 91561, 124440,
                                                                       62386, 62926, 95257,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 129592, 0, 3,
                                                                       124440, 92149, 125224,
                                                                       62926, 63466, 96013,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 130600, 0, 3,
                                                                       125224, 92737, 126008,
                                                                       63466, 64006, 96769,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 131608, 0, 3,
                                                                       126008, 93325, 126792,
                                                                       64006, 64546, 97525,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 132616, 0, 3,
                                                                       127576, 94501, 128584,
                                                                       65626, 66301, 98281,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 133876, 0, 3,
                                                                       128584, 95257, 129592,
                                                                       66301, 66976, 99226,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 135136, 0, 3,
                                                                       129592, 96013, 130600,
                                                                       66976, 67651, 100171,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 136396, 0, 3,
                                                                       130600, 96769, 131608,
                                                                       67651, 68326, 101116,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 137656, 0, 3,
                                                                       132616, 98281, 133876,
                                                                       69676, 70501, 102061,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 139196, 0, 3,
                                                                       133876, 99226, 135136,
                                                                       70501, 71326, 103216,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 140736, 0, 3,
                                                                       135136, 100171, 136396,
                                                                       71326, 72151, 104371,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 142276, 0, 3,
                                                                       137656, 102061, 139196,
                                                                       73801, 74791, 105526,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 144124, 0, 3,
                                                                       139196, 103216, 140736,
                                                                       74791, 75781, 106912,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 145972, 0, 3,
                                                                       142276, 105526, 144124,
                                                                       77761, 78931, 108298,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 148156, 122872, 784, ncols);

                    simdfunc::contract_primitives(buffer, 149304, 127576, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 150780, 132616, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 152625, 137656, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 154880, 142276, 1848, ncols);

                    simdfunc::contract_primitives(buffer, 157586, 145972, 2184, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 148940, 148156, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 150312, 149304, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 152040, 150780, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 154165, 152625, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 156728, 154880, 66, 1, nmax);

        simdtrf::transform_i_inner(buffer, 159770, 157586, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 160784, 148940, 150312, 13,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 161876, 150312, 152040, 13,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 163280, 152040, 154165, 13,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 165035, 154165, 156728, 13,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 167180, 156728, 159770, 13,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 169754, 160784, 161876, 13,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 171938, 161876, 163280, 13,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 174746, 163280, 165035, 13,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 178256, 165035, 167180, 13,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 182546, 169754, 171938, 13,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 186186, 171938, 174746, 13,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 190866, 174746, 178256, 13,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 196716, 182546, 186186, 13,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 202176, 186186, 190866, 13,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 209196, 196716, 202176, 13,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 216840, 209196, 28, 13, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 216840, 143, nmax);
    }

    for (size_t m = 0; m < 1859; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
