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


#include "SimdThreeCenterElectronRepulsionRecHHK.hpp"

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
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHG.hpp"
#include "SimdTransferHH.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hhk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hhk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 232569, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1815 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 232569, 165888, 11811, dimensions);

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
                                                        16, 17}, ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 8, 9,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 9, 10,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 91, 0, 3, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 12, 13,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 13, 14,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 109, 0, 3, 14, 15,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 115, 0, 3, 15, 16,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 121, 0, 3, 16, 17,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 127, 0, 3, 17, 18,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 18, 19,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 139, 0, 3, 19, 20,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 145, 0, 3, 20, 21,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 151, 0, 3, 21, 22,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 157, 0, 3, 22, 23,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 25, 28,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 28, 31,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 31, 34,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 34, 37,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 37, 40,
                                                                       97, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 40, 43,
                                                                       103, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 43, 46,
                                                                       109, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 46, 49,
                                                                       115, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 49, 52,
                                                                       121, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 52, 55,
                                                                       127, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 263, 0, 3, 55, 58,
                                                                       133, 139, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 58, 61,
                                                                       139, 145, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 283, 0, 3, 61, 64,
                                                                       145, 151, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 293, 0, 3, 64, 67,
                                                                       151, 157, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 303, 0, 3, 73, 79,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 79, 85,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 333, 0, 3, 85, 91,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 91, 97,
                                                                       193, 203, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 363, 0, 3, 97,
                                                                       103, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 103,
                                                                       109, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 393, 0, 3, 109,
                                                                       115, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 115,
                                                                       121, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 423, 0, 3, 121,
                                                                       127, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 127,
                                                                       133, 253, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 453, 0, 3, 133,
                                                                       139, 263, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 139,
                                                                       145, 273, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 483, 0, 3, 145,
                                                                       151, 283, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 163,
                                                                       173, 303, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 519, 0, 3, 173,
                                                                       183, 318, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 540, 0, 3, 183,
                                                                       193, 333, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 561, 0, 3, 193,
                                                                       203, 348, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 582, 0, 3, 203,
                                                                       213, 363, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 603, 0, 3, 213,
                                                                       223, 378, 393, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 624, 0, 3, 223,
                                                                       233, 393, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 645, 0, 3, 233,
                                                                       243, 408, 423, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 666, 0, 3, 243,
                                                                       253, 423, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 687, 0, 3, 253,
                                                                       263, 438, 453, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 708, 0, 3, 263,
                                                                       273, 453, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 729, 0, 3, 273,
                                                                       283, 468, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 750, 0, 3, 303,
                                                                       318, 498, 519, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 778, 0, 3, 318,
                                                                       333, 519, 540, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 806, 0, 3, 333,
                                                                       348, 540, 561, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 834, 0, 3, 348,
                                                                       363, 561, 582, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 862, 0, 3, 363,
                                                                       378, 582, 603, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 890, 0, 3, 378,
                                                                       393, 603, 624, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 918, 0, 3, 393,
                                                                       408, 624, 645, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 946, 0, 3, 408,
                                                                       423, 645, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 974, 0, 3, 423,
                                                                       438, 666, 687, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 438,
                                                                       453, 687, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 453,
                                                                       468, 708, 729, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 498,
                                                                       519, 750, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1094, 0, 3, 519,
                                                                       540, 778, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1130, 0, 3, 540,
                                                                       561, 806, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1166, 0, 3, 561,
                                                                       582, 834, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1202, 0, 3, 582,
                                                                       603, 862, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1238, 0, 3, 603,
                                                                       624, 890, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1274, 0, 3, 624,
                                                                       645, 918, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1310, 0, 3, 645,
                                                                       666, 946, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1346, 0, 3, 666,
                                                                       687, 974, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1382, 0, 3, 687,
                                                                       708, 1002, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1418, 0, 3, 750,
                                                                       778, 1058, 1094, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1463, 0, 3, 778,
                                                                       806, 1094, 1130, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1508, 0, 3, 806,
                                                                       834, 1130, 1166, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1553, 0, 3, 834,
                                                                       862, 1166, 1202, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1598, 0, 3, 862,
                                                                       890, 1202, 1238, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1643, 0, 3, 890,
                                                                       918, 1238, 1274, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 918,
                                                                       946, 1274, 1310, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1733, 0, 3, 946,
                                                                       974, 1310, 1346, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1778, 0, 3, 974,
                                                                       1002, 1346, 1382, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1823, 0, 3, 1058,
                                                                       1094, 1418, 1463, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1878, 0, 3, 1094,
                                                                       1130, 1463, 1508, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1933, 0, 3, 1130,
                                                                       1166, 1508, 1553, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1166,
                                                                       1202, 1553, 1598, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2043, 0, 3, 1202,
                                                                       1238, 1598, 1643, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2098, 0, 3, 1238,
                                                                       1274, 1643, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1274,
                                                                       1310, 1688, 1733, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2208, 0, 3, 1310,
                                                                       1346, 1733, 1778, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2263, 0, 3, 1418,
                                                                       1463, 1823, 1878, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2329, 0, 3, 1463,
                                                                       1508, 1878, 1933, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2395, 0, 3, 1508,
                                                                       1553, 1933, 1988, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2461, 0, 3, 1553,
                                                                       1598, 1988, 2043, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2527, 0, 3, 1598,
                                                                       1643, 2043, 2098, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2593, 0, 3, 1643,
                                                                       1688, 2098, 2153, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2659, 0, 3, 1688,
                                                                       1733, 2153, 2208, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2725, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2728, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2731, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2734, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2737, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2740, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2743, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2746, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2749, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2752, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2755, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2758, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2761, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2764, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2767, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2770, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2773, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2776, 3, 8, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2785, 3, 9, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2794, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2803, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2812, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2821, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2830, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2839, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2848, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2857, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2866, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2875, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2884, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2893, 3, 21, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2902, 3, 22, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2911, 3, 23, 70,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2920, 3, 25, 73,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2938, 3, 28, 79,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2956, 3, 31, 85,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2974, 3, 34, 91,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2992, 3, 37, 97,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3010, 3, 40, 103,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3028, 3, 43, 109,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3046, 3, 46, 115,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3064, 3, 49, 121,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3082, 3, 52, 127,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3100, 3, 55, 133,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3118, 3, 58, 139,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3136, 3, 61, 145,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3154, 3, 64, 151,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3172, 3, 67, 157,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3190, 3, 73, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3220, 3, 79, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3250, 3, 85, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3280, 3, 91, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3310, 3, 97, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3340, 3, 103, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3370, 3, 109, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3400, 3, 115, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3430, 3, 121, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3460, 3, 127, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3490, 3, 133, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3520, 3, 139, 273,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3550, 3, 145, 283,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3580, 3, 151, 293,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3610, 3, 163, 303,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3655, 3, 173, 318,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3700, 3, 183, 333,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3745, 3, 193, 348,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3790, 3, 203, 363,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3835, 3, 213, 378,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3880, 3, 223, 393,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3925, 3, 233, 408,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3970, 3, 243, 423,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4015, 3, 253, 438,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4060, 3, 263, 453,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4105, 3, 273, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4150, 3, 283, 483,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4195, 3, 303, 498,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4258, 3, 318, 519,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4321, 3, 333, 540,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4384, 3, 348, 561,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4447, 3, 363, 582,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4510, 3, 378, 603,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4573, 3, 393, 624,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4636, 3, 408, 645,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4699, 3, 423, 666,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4762, 3, 438, 687,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4825, 3, 453, 708,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4888, 3, 468, 729,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4951, 3, 498, 750,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5035, 3, 519, 778,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5119, 3, 540, 806,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5203, 3, 561, 834,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5287, 3, 582, 862,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5371, 3, 603, 890,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5455, 3, 624, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5539, 3, 645, 946,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5623, 3, 666, 974,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5707, 3, 687,
                                                                       1002, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5791, 3, 708,
                                                                       1030, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5875, 3, 750,
                                                                       1058, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5983, 3, 778,
                                                                       1094, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6091, 3, 806,
                                                                       1130, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6199, 3, 834,
                                                                       1166, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6307, 3, 862,
                                                                       1202, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6415, 3, 890,
                                                                       1238, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6523, 3, 918,
                                                                       1274, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6631, 3, 946,
                                                                       1310, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6739, 3, 974,
                                                                       1346, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6847, 3, 1002,
                                                                       1382, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6955, 3, 1058,
                                                                       1418, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7090, 3, 1094,
                                                                       1463, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7225, 3, 1130,
                                                                       1508, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7360, 3, 1166,
                                                                       1553, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7495, 3, 1202,
                                                                       1598, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7630, 3, 1238,
                                                                       1643, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7765, 3, 1274,
                                                                       1688, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7900, 3, 1310,
                                                                       1733, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8035, 3, 1346,
                                                                       1778, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8170, 3, 1418,
                                                                       1823, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8335, 3, 1463,
                                                                       1878, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8500, 3, 1508,
                                                                       1933, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8665, 3, 1553,
                                                                       1988, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8830, 3, 1598,
                                                                       2043, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8995, 3, 1643,
                                                                       2098, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9160, 3, 1688,
                                                                       2153, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9325, 3, 1733,
                                                                       2208, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9490, 3, 1823,
                                                                       2263, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9688, 3, 1878,
                                                                       2329, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9886, 3, 1933,
                                                                       2395, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10084, 3, 1988,
                                                                       2461, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10282, 3, 2043,
                                                                       2527, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10480, 3, 2098,
                                                                       2593, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10678, 3, 2153,
                                                                       2659, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10876, 3, 8, 9,
                                                                       2731, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10882, 3, 9, 10,
                                                                       2734, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10888, 3, 10, 11,
                                                                       2737, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10894, 3, 11, 12,
                                                                       2740, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10900, 3, 12, 13,
                                                                       2743, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10906, 3, 13, 14,
                                                                       2746, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10912, 3, 14, 15,
                                                                       2749, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10918, 3, 15, 16,
                                                                       2752, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10924, 3, 16, 17,
                                                                       2755, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10930, 3, 17, 18,
                                                                       2758, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10936, 3, 18, 19,
                                                                       2761, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10942, 3, 19, 20,
                                                                       2764, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10948, 3, 20, 21,
                                                                       2767, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10954, 3, 21, 22,
                                                                       2770, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10960, 3, 22, 23,
                                                                       2773, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10966, 0, 3,
                                                                       10876, 2731, 10882, 2794,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10984, 0, 3,
                                                                       10882, 2734, 10888, 2803,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11002, 0, 3,
                                                                       10888, 2737, 10894, 2812,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11020, 0, 3,
                                                                       10894, 2740, 10900, 2821,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11038, 0, 3,
                                                                       10900, 2743, 10906, 2830,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11056, 0, 3,
                                                                       10906, 2746, 10912, 2839,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11074, 0, 3,
                                                                       10912, 2749, 10918, 2848,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11092, 0, 3,
                                                                       10918, 2752, 10924, 2857,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11110, 0, 3,
                                                                       10924, 2755, 10930, 2866,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11128, 0, 3,
                                                                       10930, 2758, 10936, 2875,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11146, 0, 3,
                                                                       10936, 2761, 10942, 2884,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11164, 0, 3,
                                                                       10942, 2764, 10948, 2893,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11182, 0, 3,
                                                                       10948, 2767, 10954, 2902,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11200, 0, 3,
                                                                       10954, 2770, 10960, 2911,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11218, 0, 3,
                                                                       10966, 2794, 10984, 73,
                                                                       79, 2956, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11254, 0, 3,
                                                                       10984, 2803, 11002, 79,
                                                                       85, 2974, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11290, 0, 3,
                                                                       11002, 2812, 11020, 85,
                                                                       91, 2992, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11326, 0, 3,
                                                                       11020, 2821, 11038, 91,
                                                                       97, 3010, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11362, 0, 3,
                                                                       11038, 2830, 11056, 97,
                                                                       103, 3028, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11398, 0, 3,
                                                                       11056, 2839, 11074, 103,
                                                                       109, 3046, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11434, 0, 3,
                                                                       11074, 2848, 11092, 109,
                                                                       115, 3064, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11470, 0, 3,
                                                                       11092, 2857, 11110, 115,
                                                                       121, 3082, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11506, 0, 3,
                                                                       11110, 2866, 11128, 121,
                                                                       127, 3100, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11542, 0, 3,
                                                                       11128, 2875, 11146, 127,
                                                                       133, 3118, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11578, 0, 3,
                                                                       11146, 2884, 11164, 133,
                                                                       139, 3136, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11614, 0, 3,
                                                                       11164, 2893, 11182, 139,
                                                                       145, 3154, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11650, 0, 3,
                                                                       11182, 2902, 11200, 145,
                                                                       151, 3172, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11686, 0, 3,
                                                                       11218, 2956, 11254, 163,
                                                                       173, 3250, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11746, 0, 3,
                                                                       11254, 2974, 11290, 173,
                                                                       183, 3280, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11806, 0, 3,
                                                                       11290, 2992, 11326, 183,
                                                                       193, 3310, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11866, 0, 3,
                                                                       11326, 3010, 11362, 193,
                                                                       203, 3340, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11926, 0, 3,
                                                                       11362, 3028, 11398, 203,
                                                                       213, 3370, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11986, 0, 3,
                                                                       11398, 3046, 11434, 213,
                                                                       223, 3400, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12046, 0, 3,
                                                                       11434, 3064, 11470, 223,
                                                                       233, 3430, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12106, 0, 3,
                                                                       11470, 3082, 11506, 233,
                                                                       243, 3460, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12166, 0, 3,
                                                                       11506, 3100, 11542, 243,
                                                                       253, 3490, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12226, 0, 3,
                                                                       11542, 3118, 11578, 253,
                                                                       263, 3520, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12286, 0, 3,
                                                                       11578, 3136, 11614, 263,
                                                                       273, 3550, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12346, 0, 3,
                                                                       11614, 3154, 11650, 273,
                                                                       283, 3580, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12406, 0, 3,
                                                                       11686, 3250, 11746, 303,
                                                                       318, 3700, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12496, 0, 3,
                                                                       11746, 3280, 11806, 318,
                                                                       333, 3745, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12586, 0, 3,
                                                                       11806, 3310, 11866, 333,
                                                                       348, 3790, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12676, 0, 3,
                                                                       11866, 3340, 11926, 348,
                                                                       363, 3835, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12766, 0, 3,
                                                                       11926, 3370, 11986, 363,
                                                                       378, 3880, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12856, 0, 3,
                                                                       11986, 3400, 12046, 378,
                                                                       393, 3925, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12946, 0, 3,
                                                                       12046, 3430, 12106, 393,
                                                                       408, 3970, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13036, 0, 3,
                                                                       12106, 3460, 12166, 408,
                                                                       423, 4015, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13126, 0, 3,
                                                                       12166, 3490, 12226, 423,
                                                                       438, 4060, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13216, 0, 3,
                                                                       12226, 3520, 12286, 438,
                                                                       453, 4105, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13306, 0, 3,
                                                                       12286, 3550, 12346, 453,
                                                                       468, 4150, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13396, 0, 3,
                                                                       12406, 3700, 12496, 498,
                                                                       519, 4321, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13522, 0, 3,
                                                                       12496, 3745, 12586, 519,
                                                                       540, 4384, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13648, 0, 3,
                                                                       12586, 3790, 12676, 540,
                                                                       561, 4447, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13774, 0, 3,
                                                                       12676, 3835, 12766, 561,
                                                                       582, 4510, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13900, 0, 3,
                                                                       12766, 3880, 12856, 582,
                                                                       603, 4573, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14026, 0, 3,
                                                                       12856, 3925, 12946, 603,
                                                                       624, 4636, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14152, 0, 3,
                                                                       12946, 3970, 13036, 624,
                                                                       645, 4699, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14278, 0, 3,
                                                                       13036, 4015, 13126, 645,
                                                                       666, 4762, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14404, 0, 3,
                                                                       13126, 4060, 13216, 666,
                                                                       687, 4825, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14530, 0, 3,
                                                                       13216, 4105, 13306, 687,
                                                                       708, 4888, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14656, 0, 3,
                                                                       13396, 4321, 13522, 750,
                                                                       778, 5119, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14824, 0, 3,
                                                                       13522, 4384, 13648, 778,
                                                                       806, 5203, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14992, 0, 3,
                                                                       13648, 4447, 13774, 806,
                                                                       834, 5287, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15160, 0, 3,
                                                                       13774, 4510, 13900, 834,
                                                                       862, 5371, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15328, 0, 3,
                                                                       13900, 4573, 14026, 862,
                                                                       890, 5455, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15496, 0, 3,
                                                                       14026, 4636, 14152, 890,
                                                                       918, 5539, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15664, 0, 3,
                                                                       14152, 4699, 14278, 918,
                                                                       946, 5623, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15832, 0, 3,
                                                                       14278, 4762, 14404, 946,
                                                                       974, 5707, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16000, 0, 3,
                                                                       14404, 4825, 14530, 974,
                                                                       1002, 5791, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16168, 0, 3,
                                                                       14656, 5119, 14824, 1058,
                                                                       1094, 6091, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16384, 0, 3,
                                                                       14824, 5203, 14992, 1094,
                                                                       1130, 6199, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16600, 0, 3,
                                                                       14992, 5287, 15160, 1130,
                                                                       1166, 6307, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16816, 0, 3,
                                                                       15160, 5371, 15328, 1166,
                                                                       1202, 6415, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17032, 0, 3,
                                                                       15328, 5455, 15496, 1202,
                                                                       1238, 6523, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17248, 0, 3,
                                                                       15496, 5539, 15664, 1238,
                                                                       1274, 6631, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17464, 0, 3,
                                                                       15664, 5623, 15832, 1274,
                                                                       1310, 6739, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17680, 0, 3,
                                                                       15832, 5707, 16000, 1310,
                                                                       1346, 6847, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17896, 0, 3,
                                                                       16168, 6091, 16384, 1418,
                                                                       1463, 7225, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18166, 0, 3,
                                                                       16384, 6199, 16600, 1463,
                                                                       1508, 7360, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18436, 0, 3,
                                                                       16600, 6307, 16816, 1508,
                                                                       1553, 7495, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18706, 0, 3,
                                                                       16816, 6415, 17032, 1553,
                                                                       1598, 7630, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18976, 0, 3,
                                                                       17032, 6523, 17248, 1598,
                                                                       1643, 7765, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19246, 0, 3,
                                                                       17248, 6631, 17464, 1643,
                                                                       1688, 7900, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19516, 0, 3,
                                                                       17464, 6739, 17680, 1688,
                                                                       1733, 8035, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19786, 0, 3,
                                                                       17896, 7225, 18166, 1823,
                                                                       1878, 8500, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20116, 0, 3,
                                                                       18166, 7360, 18436, 1878,
                                                                       1933, 8665, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20446, 0, 3,
                                                                       18436, 7495, 18706, 1933,
                                                                       1988, 8830, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20776, 0, 3,
                                                                       18706, 7630, 18976, 1988,
                                                                       2043, 8995, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 21106, 0, 3,
                                                                       18976, 7765, 19246, 2043,
                                                                       2098, 9160, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 21436, 0, 3,
                                                                       19246, 7900, 19516, 2098,
                                                                       2153, 9325, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 21766, 0, 3,
                                                                       19786, 8500, 20116, 2263,
                                                                       2329, 9886, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 22162, 0, 3,
                                                                       20116, 8665, 20446, 2329,
                                                                       2395, 10084, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 22558, 0, 3,
                                                                       20446, 8830, 20776, 2395,
                                                                       2461, 10282, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 22954, 0, 3,
                                                                       20776, 8995, 21106, 2461,
                                                                       2527, 10480, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 23350, 0, 3,
                                                                       21106, 9160, 21436, 2527,
                                                                       2593, 10678, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23746, 3, 2725,
                                                                       2728, 10876, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23756, 3, 2728,
                                                                       2731, 10882, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23766, 3, 2731,
                                                                       2734, 10888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23776, 3, 2734,
                                                                       2737, 10894, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23786, 3, 2737,
                                                                       2740, 10900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23796, 3, 2740,
                                                                       2743, 10906, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23806, 3, 2743,
                                                                       2746, 10912, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23816, 3, 2746,
                                                                       2749, 10918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23826, 3, 2749,
                                                                       2752, 10924, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23836, 3, 2752,
                                                                       2755, 10930, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23846, 3, 2755,
                                                                       2758, 10936, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23856, 3, 2758,
                                                                       2761, 10942, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23866, 3, 2761,
                                                                       2764, 10948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23876, 3, 2764,
                                                                       2767, 10954, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23886, 3, 2767,
                                                                       2770, 10960, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 23896, 0, 3,
                                                                       23746, 10876, 23756, 2776,
                                                                       2785, 10966, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 23926, 0, 3,
                                                                       23756, 10882, 23766, 2785,
                                                                       2794, 10984, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 23956, 0, 3,
                                                                       23766, 10888, 23776, 2794,
                                                                       2803, 11002, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 23986, 0, 3,
                                                                       23776, 10894, 23786, 2803,
                                                                       2812, 11020, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24016, 0, 3,
                                                                       23786, 10900, 23796, 2812,
                                                                       2821, 11038, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24046, 0, 3,
                                                                       23796, 10906, 23806, 2821,
                                                                       2830, 11056, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24076, 0, 3,
                                                                       23806, 10912, 23816, 2830,
                                                                       2839, 11074, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24106, 0, 3,
                                                                       23816, 10918, 23826, 2839,
                                                                       2848, 11092, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24136, 0, 3,
                                                                       23826, 10924, 23836, 2848,
                                                                       2857, 11110, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24166, 0, 3,
                                                                       23836, 10930, 23846, 2857,
                                                                       2866, 11128, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24196, 0, 3,
                                                                       23846, 10936, 23856, 2866,
                                                                       2875, 11146, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24226, 0, 3,
                                                                       23856, 10942, 23866, 2875,
                                                                       2884, 11164, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24256, 0, 3,
                                                                       23866, 10948, 23876, 2884,
                                                                       2893, 11182, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24286, 0, 3,
                                                                       23876, 10954, 23886, 2893,
                                                                       2902, 11200, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24316, 0, 3,
                                                                       23896, 10966, 23926, 2920,
                                                                       2938, 11218, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24376, 0, 3,
                                                                       23926, 10984, 23956, 2938,
                                                                       2956, 11254, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24436, 0, 3,
                                                                       23956, 11002, 23986, 2956,
                                                                       2974, 11290, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24496, 0, 3,
                                                                       23986, 11020, 24016, 2974,
                                                                       2992, 11326, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24556, 0, 3,
                                                                       24016, 11038, 24046, 2992,
                                                                       3010, 11362, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24616, 0, 3,
                                                                       24046, 11056, 24076, 3010,
                                                                       3028, 11398, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24676, 0, 3,
                                                                       24076, 11074, 24106, 3028,
                                                                       3046, 11434, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24736, 0, 3,
                                                                       24106, 11092, 24136, 3046,
                                                                       3064, 11470, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24796, 0, 3,
                                                                       24136, 11110, 24166, 3064,
                                                                       3082, 11506, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24856, 0, 3,
                                                                       24166, 11128, 24196, 3082,
                                                                       3100, 11542, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24916, 0, 3,
                                                                       24196, 11146, 24226, 3100,
                                                                       3118, 11578, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24976, 0, 3,
                                                                       24226, 11164, 24256, 3118,
                                                                       3136, 11614, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25036, 0, 3,
                                                                       24256, 11182, 24286, 3136,
                                                                       3154, 11650, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25096, 0, 3,
                                                                       24316, 11218, 24376, 3190,
                                                                       3220, 11686, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25196, 0, 3,
                                                                       24376, 11254, 24436, 3220,
                                                                       3250, 11746, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25296, 0, 3,
                                                                       24436, 11290, 24496, 3250,
                                                                       3280, 11806, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25396, 0, 3,
                                                                       24496, 11326, 24556, 3280,
                                                                       3310, 11866, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25496, 0, 3,
                                                                       24556, 11362, 24616, 3310,
                                                                       3340, 11926, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25596, 0, 3,
                                                                       24616, 11398, 24676, 3340,
                                                                       3370, 11986, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25696, 0, 3,
                                                                       24676, 11434, 24736, 3370,
                                                                       3400, 12046, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25796, 0, 3,
                                                                       24736, 11470, 24796, 3400,
                                                                       3430, 12106, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25896, 0, 3,
                                                                       24796, 11506, 24856, 3430,
                                                                       3460, 12166, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25996, 0, 3,
                                                                       24856, 11542, 24916, 3460,
                                                                       3490, 12226, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26096, 0, 3,
                                                                       24916, 11578, 24976, 3490,
                                                                       3520, 12286, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26196, 0, 3,
                                                                       24976, 11614, 25036, 3520,
                                                                       3550, 12346, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26296, 0, 3,
                                                                       25096, 11686, 25196, 3610,
                                                                       3655, 12406, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26446, 0, 3,
                                                                       25196, 11746, 25296, 3655,
                                                                       3700, 12496, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26596, 0, 3,
                                                                       25296, 11806, 25396, 3700,
                                                                       3745, 12586, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26746, 0, 3,
                                                                       25396, 11866, 25496, 3745,
                                                                       3790, 12676, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26896, 0, 3,
                                                                       25496, 11926, 25596, 3790,
                                                                       3835, 12766, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27046, 0, 3,
                                                                       25596, 11986, 25696, 3835,
                                                                       3880, 12856, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27196, 0, 3,
                                                                       25696, 12046, 25796, 3880,
                                                                       3925, 12946, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27346, 0, 3,
                                                                       25796, 12106, 25896, 3925,
                                                                       3970, 13036, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27496, 0, 3,
                                                                       25896, 12166, 25996, 3970,
                                                                       4015, 13126, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27646, 0, 3,
                                                                       25996, 12226, 26096, 4015,
                                                                       4060, 13216, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27796, 0, 3,
                                                                       26096, 12286, 26196, 4060,
                                                                       4105, 13306, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 27946, 0, 3,
                                                                       26296, 12406, 26446, 4195,
                                                                       4258, 13396, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28156, 0, 3,
                                                                       26446, 12496, 26596, 4258,
                                                                       4321, 13522, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28366, 0, 3,
                                                                       26596, 12586, 26746, 4321,
                                                                       4384, 13648, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28576, 0, 3,
                                                                       26746, 12676, 26896, 4384,
                                                                       4447, 13774, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28786, 0, 3,
                                                                       26896, 12766, 27046, 4447,
                                                                       4510, 13900, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28996, 0, 3,
                                                                       27046, 12856, 27196, 4510,
                                                                       4573, 14026, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29206, 0, 3,
                                                                       27196, 12946, 27346, 4573,
                                                                       4636, 14152, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29416, 0, 3,
                                                                       27346, 13036, 27496, 4636,
                                                                       4699, 14278, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29626, 0, 3,
                                                                       27496, 13126, 27646, 4699,
                                                                       4762, 14404, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29836, 0, 3,
                                                                       27646, 13216, 27796, 4762,
                                                                       4825, 14530, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30046, 0, 3,
                                                                       27946, 13396, 28156, 4951,
                                                                       5035, 14656, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30326, 0, 3,
                                                                       28156, 13522, 28366, 5035,
                                                                       5119, 14824, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30606, 0, 3,
                                                                       28366, 13648, 28576, 5119,
                                                                       5203, 14992, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30886, 0, 3,
                                                                       28576, 13774, 28786, 5203,
                                                                       5287, 15160, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 31166, 0, 3,
                                                                       28786, 13900, 28996, 5287,
                                                                       5371, 15328, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 31446, 0, 3,
                                                                       28996, 14026, 29206, 5371,
                                                                       5455, 15496, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 31726, 0, 3,
                                                                       29206, 14152, 29416, 5455,
                                                                       5539, 15664, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 32006, 0, 3,
                                                                       29416, 14278, 29626, 5539,
                                                                       5623, 15832, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 32286, 0, 3,
                                                                       29626, 14404, 29836, 5623,
                                                                       5707, 16000, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32566, 0, 3,
                                                                       30046, 14656, 30326, 5875,
                                                                       5983, 16168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32926, 0, 3,
                                                                       30326, 14824, 30606, 5983,
                                                                       6091, 16384, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33286, 0, 3,
                                                                       30606, 14992, 30886, 6091,
                                                                       6199, 16600, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33646, 0, 3,
                                                                       30886, 15160, 31166, 6199,
                                                                       6307, 16816, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 34006, 0, 3,
                                                                       31166, 15328, 31446, 6307,
                                                                       6415, 17032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 34366, 0, 3,
                                                                       31446, 15496, 31726, 6415,
                                                                       6523, 17248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 34726, 0, 3,
                                                                       31726, 15664, 32006, 6523,
                                                                       6631, 17464, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 35086, 0, 3,
                                                                       32006, 15832, 32286, 6631,
                                                                       6739, 17680, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 35446, 0, 3,
                                                                       32566, 16168, 32926, 6955,
                                                                       7090, 17896, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 35896, 0, 3,
                                                                       32926, 16384, 33286, 7090,
                                                                       7225, 18166, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 36346, 0, 3,
                                                                       33286, 16600, 33646, 7225,
                                                                       7360, 18436, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 36796, 0, 3,
                                                                       33646, 16816, 34006, 7360,
                                                                       7495, 18706, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 37246, 0, 3,
                                                                       34006, 17032, 34366, 7495,
                                                                       7630, 18976, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 37696, 0, 3,
                                                                       34366, 17248, 34726, 7630,
                                                                       7765, 19246, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 38146, 0, 3,
                                                                       34726, 17464, 35086, 7765,
                                                                       7900, 19516, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 38596, 0, 3,
                                                                       35446, 17896, 35896, 8170,
                                                                       8335, 19786, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 39146, 0, 3,
                                                                       35896, 18166, 36346, 8335,
                                                                       8500, 20116, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 39696, 0, 3,
                                                                       36346, 18436, 36796, 8500,
                                                                       8665, 20446, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 40246, 0, 3,
                                                                       36796, 18706, 37246, 8665,
                                                                       8830, 20776, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 40796, 0, 3,
                                                                       37246, 18976, 37696, 8830,
                                                                       8995, 21106, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 41346, 0, 3,
                                                                       37696, 19246, 38146, 8995,
                                                                       9160, 21436, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 41896, 0, 3,
                                                                       38596, 19786, 39146, 9490,
                                                                       9688, 21766, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 42556, 0, 3,
                                                                       39146, 20116, 39696, 9688,
                                                                       9886, 22162, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 43216, 0, 3,
                                                                       39696, 20446, 40246, 9886,
                                                                       10084, 22558, ncols,
                                                                       gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 43876, 0, 3,
                                                                       40246, 20776, 40796,
                                                                       10084, 10282, 22954,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 44536, 0, 3,
                                                                       40796, 21106, 41346,
                                                                       10282, 10480, 23350,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45196, 3, 10876,
                                                                       10882, 23766, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45211, 3, 10882,
                                                                       10888, 23776, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45226, 3, 10888,
                                                                       10894, 23786, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45241, 3, 10894,
                                                                       10900, 23796, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45256, 3, 10900,
                                                                       10906, 23806, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45271, 3, 10906,
                                                                       10912, 23816, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45286, 3, 10912,
                                                                       10918, 23826, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45301, 3, 10918,
                                                                       10924, 23836, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45316, 3, 10924,
                                                                       10930, 23846, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45331, 3, 10930,
                                                                       10936, 23856, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45346, 3, 10936,
                                                                       10942, 23866, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45361, 3, 10942,
                                                                       10948, 23876, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45376, 3, 10948,
                                                                       10954, 23886, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45391, 0, 3,
                                                                       45196, 23766, 45211,
                                                                       10966, 10984, 23956,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45436, 0, 3,
                                                                       45211, 23776, 45226,
                                                                       10984, 11002, 23986,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45481, 0, 3,
                                                                       45226, 23786, 45241,
                                                                       11002, 11020, 24016,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45526, 0, 3,
                                                                       45241, 23796, 45256,
                                                                       11020, 11038, 24046,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45571, 0, 3,
                                                                       45256, 23806, 45271,
                                                                       11038, 11056, 24076,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45616, 0, 3,
                                                                       45271, 23816, 45286,
                                                                       11056, 11074, 24106,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45661, 0, 3,
                                                                       45286, 23826, 45301,
                                                                       11074, 11092, 24136,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45706, 0, 3,
                                                                       45301, 23836, 45316,
                                                                       11092, 11110, 24166,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45751, 0, 3,
                                                                       45316, 23846, 45331,
                                                                       11110, 11128, 24196,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45796, 0, 3,
                                                                       45331, 23856, 45346,
                                                                       11128, 11146, 24226,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45841, 0, 3,
                                                                       45346, 23866, 45361,
                                                                       11146, 11164, 24256,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45886, 0, 3,
                                                                       45361, 23876, 45376,
                                                                       11164, 11182, 24286,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45931, 0, 3,
                                                                       45391, 23956, 45436,
                                                                       11218, 11254, 24436,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46021, 0, 3,
                                                                       45436, 23986, 45481,
                                                                       11254, 11290, 24496,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46111, 0, 3,
                                                                       45481, 24016, 45526,
                                                                       11290, 11326, 24556,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46201, 0, 3,
                                                                       45526, 24046, 45571,
                                                                       11326, 11362, 24616,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46291, 0, 3,
                                                                       45571, 24076, 45616,
                                                                       11362, 11398, 24676,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46381, 0, 3,
                                                                       45616, 24106, 45661,
                                                                       11398, 11434, 24736,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46471, 0, 3,
                                                                       45661, 24136, 45706,
                                                                       11434, 11470, 24796,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46561, 0, 3,
                                                                       45706, 24166, 45751,
                                                                       11470, 11506, 24856,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46651, 0, 3,
                                                                       45751, 24196, 45796,
                                                                       11506, 11542, 24916,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46741, 0, 3,
                                                                       45796, 24226, 45841,
                                                                       11542, 11578, 24976,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46831, 0, 3,
                                                                       45841, 24256, 45886,
                                                                       11578, 11614, 25036,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46921, 0, 3,
                                                                       45931, 24436, 46021,
                                                                       11686, 11746, 25296,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47071, 0, 3,
                                                                       46021, 24496, 46111,
                                                                       11746, 11806, 25396,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47221, 0, 3,
                                                                       46111, 24556, 46201,
                                                                       11806, 11866, 25496,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47371, 0, 3,
                                                                       46201, 24616, 46291,
                                                                       11866, 11926, 25596,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47521, 0, 3,
                                                                       46291, 24676, 46381,
                                                                       11926, 11986, 25696,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47671, 0, 3,
                                                                       46381, 24736, 46471,
                                                                       11986, 12046, 25796,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47821, 0, 3,
                                                                       46471, 24796, 46561,
                                                                       12046, 12106, 25896,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47971, 0, 3,
                                                                       46561, 24856, 46651,
                                                                       12106, 12166, 25996,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48121, 0, 3,
                                                                       46651, 24916, 46741,
                                                                       12166, 12226, 26096,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48271, 0, 3,
                                                                       46741, 24976, 46831,
                                                                       12226, 12286, 26196,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48421, 0, 3,
                                                                       46921, 25296, 47071,
                                                                       12406, 12496, 26596,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48646, 0, 3,
                                                                       47071, 25396, 47221,
                                                                       12496, 12586, 26746,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48871, 0, 3,
                                                                       47221, 25496, 47371,
                                                                       12586, 12676, 26896,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49096, 0, 3,
                                                                       47371, 25596, 47521,
                                                                       12676, 12766, 27046,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49321, 0, 3,
                                                                       47521, 25696, 47671,
                                                                       12766, 12856, 27196,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49546, 0, 3,
                                                                       47671, 25796, 47821,
                                                                       12856, 12946, 27346,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49771, 0, 3,
                                                                       47821, 25896, 47971,
                                                                       12946, 13036, 27496,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49996, 0, 3,
                                                                       47971, 25996, 48121,
                                                                       13036, 13126, 27646,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50221, 0, 3,
                                                                       48121, 26096, 48271,
                                                                       13126, 13216, 27796,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50446, 0, 3,
                                                                       48421, 26596, 48646,
                                                                       13396, 13522, 28366,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50761, 0, 3,
                                                                       48646, 26746, 48871,
                                                                       13522, 13648, 28576,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51076, 0, 3,
                                                                       48871, 26896, 49096,
                                                                       13648, 13774, 28786,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51391, 0, 3,
                                                                       49096, 27046, 49321,
                                                                       13774, 13900, 28996,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51706, 0, 3,
                                                                       49321, 27196, 49546,
                                                                       13900, 14026, 29206,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52021, 0, 3,
                                                                       49546, 27346, 49771,
                                                                       14026, 14152, 29416,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52336, 0, 3,
                                                                       49771, 27496, 49996,
                                                                       14152, 14278, 29626,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52651, 0, 3,
                                                                       49996, 27646, 50221,
                                                                       14278, 14404, 29836,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 52966, 0, 3,
                                                                       50446, 28366, 50761,
                                                                       14656, 14824, 30606,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 53386, 0, 3,
                                                                       50761, 28576, 51076,
                                                                       14824, 14992, 30886,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 53806, 0, 3,
                                                                       51076, 28786, 51391,
                                                                       14992, 15160, 31166,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 54226, 0, 3,
                                                                       51391, 28996, 51706,
                                                                       15160, 15328, 31446,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 54646, 0, 3,
                                                                       51706, 29206, 52021,
                                                                       15328, 15496, 31726,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55066, 0, 3,
                                                                       52021, 29416, 52336,
                                                                       15496, 15664, 32006,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55486, 0, 3,
                                                                       52336, 29626, 52651,
                                                                       15664, 15832, 32286,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 55906, 0, 3,
                                                                       52966, 30606, 53386,
                                                                       16168, 16384, 33286,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 56446, 0, 3,
                                                                       53386, 30886, 53806,
                                                                       16384, 16600, 33646,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 56986, 0, 3,
                                                                       53806, 31166, 54226,
                                                                       16600, 16816, 34006,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 57526, 0, 3,
                                                                       54226, 31446, 54646,
                                                                       16816, 17032, 34366,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 58066, 0, 3,
                                                                       54646, 31726, 55066,
                                                                       17032, 17248, 34726,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 58606, 0, 3,
                                                                       55066, 32006, 55486,
                                                                       17248, 17464, 35086,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 59146, 0, 3,
                                                                       55906, 33286, 56446,
                                                                       17896, 18166, 36346,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 59821, 0, 3,
                                                                       56446, 33646, 56986,
                                                                       18166, 18436, 36796,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 60496, 0, 3,
                                                                       56986, 34006, 57526,
                                                                       18436, 18706, 37246,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 61171, 0, 3,
                                                                       57526, 34366, 58066,
                                                                       18706, 18976, 37696,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 61846, 0, 3,
                                                                       58066, 34726, 58606,
                                                                       18976, 19246, 38146,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 62521, 0, 3,
                                                                       59146, 36346, 59821,
                                                                       19786, 20116, 39696,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 63346, 0, 3,
                                                                       59821, 36796, 60496,
                                                                       20116, 20446, 40246,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 64171, 0, 3,
                                                                       60496, 37246, 61171,
                                                                       20446, 20776, 40796,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 64996, 0, 3,
                                                                       61171, 37696, 61846,
                                                                       20776, 21106, 41346,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 65821, 0, 3,
                                                                       62521, 39696, 63346,
                                                                       21766, 22162, 43216,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 66811, 0, 3,
                                                                       63346, 40246, 64171,
                                                                       22162, 22558, 43876,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 67801, 0, 3,
                                                                       64171, 40796, 64996,
                                                                       22558, 22954, 44536,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68791, 3, 23746,
                                                                       23756, 45196, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68812, 3, 23756,
                                                                       23766, 45211, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68833, 3, 23766,
                                                                       23776, 45226, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68854, 3, 23776,
                                                                       23786, 45241, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68875, 3, 23786,
                                                                       23796, 45256, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68896, 3, 23796,
                                                                       23806, 45271, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68917, 3, 23806,
                                                                       23816, 45286, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68938, 3, 23816,
                                                                       23826, 45301, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68959, 3, 23826,
                                                                       23836, 45316, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68980, 3, 23836,
                                                                       23846, 45331, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 69001, 3, 23846,
                                                                       23856, 45346, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 69022, 3, 23856,
                                                                       23866, 45361, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 69043, 3, 23866,
                                                                       23876, 45376, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69064, 0, 3,
                                                                       68791, 45196, 68812,
                                                                       23896, 23926, 45391,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69127, 0, 3,
                                                                       68812, 45211, 68833,
                                                                       23926, 23956, 45436,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69190, 0, 3,
                                                                       68833, 45226, 68854,
                                                                       23956, 23986, 45481,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69253, 0, 3,
                                                                       68854, 45241, 68875,
                                                                       23986, 24016, 45526,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69316, 0, 3,
                                                                       68875, 45256, 68896,
                                                                       24016, 24046, 45571,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69379, 0, 3,
                                                                       68896, 45271, 68917,
                                                                       24046, 24076, 45616,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69442, 0, 3,
                                                                       68917, 45286, 68938,
                                                                       24076, 24106, 45661,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69505, 0, 3,
                                                                       68938, 45301, 68959,
                                                                       24106, 24136, 45706,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69568, 0, 3,
                                                                       68959, 45316, 68980,
                                                                       24136, 24166, 45751,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69631, 0, 3,
                                                                       68980, 45331, 69001,
                                                                       24166, 24196, 45796,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69694, 0, 3,
                                                                       69001, 45346, 69022,
                                                                       24196, 24226, 45841,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69757, 0, 3,
                                                                       69022, 45361, 69043,
                                                                       24226, 24256, 45886,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 69820, 0, 3,
                                                                       69064, 45391, 69127,
                                                                       24316, 24376, 45931,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 69946, 0, 3,
                                                                       69127, 45436, 69190,
                                                                       24376, 24436, 46021,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70072, 0, 3,
                                                                       69190, 45481, 69253,
                                                                       24436, 24496, 46111,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70198, 0, 3,
                                                                       69253, 45526, 69316,
                                                                       24496, 24556, 46201,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70324, 0, 3,
                                                                       69316, 45571, 69379,
                                                                       24556, 24616, 46291,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70450, 0, 3,
                                                                       69379, 45616, 69442,
                                                                       24616, 24676, 46381,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70576, 0, 3,
                                                                       69442, 45661, 69505,
                                                                       24676, 24736, 46471,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70702, 0, 3,
                                                                       69505, 45706, 69568,
                                                                       24736, 24796, 46561,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70828, 0, 3,
                                                                       69568, 45751, 69631,
                                                                       24796, 24856, 46651,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70954, 0, 3,
                                                                       69631, 45796, 69694,
                                                                       24856, 24916, 46741,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 71080, 0, 3,
                                                                       69694, 45841, 69757,
                                                                       24916, 24976, 46831,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 71206, 0, 3,
                                                                       69820, 45931, 69946,
                                                                       25096, 25196, 46921,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 71416, 0, 3,
                                                                       69946, 46021, 70072,
                                                                       25196, 25296, 47071,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 71626, 0, 3,
                                                                       70072, 46111, 70198,
                                                                       25296, 25396, 47221,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 71836, 0, 3,
                                                                       70198, 46201, 70324,
                                                                       25396, 25496, 47371,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 72046, 0, 3,
                                                                       70324, 46291, 70450,
                                                                       25496, 25596, 47521,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 72256, 0, 3,
                                                                       70450, 46381, 70576,
                                                                       25596, 25696, 47671,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 72466, 0, 3,
                                                                       70576, 46471, 70702,
                                                                       25696, 25796, 47821,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 72676, 0, 3,
                                                                       70702, 46561, 70828,
                                                                       25796, 25896, 47971,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 72886, 0, 3,
                                                                       70828, 46651, 70954,
                                                                       25896, 25996, 48121,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 73096, 0, 3,
                                                                       70954, 46741, 71080,
                                                                       25996, 26096, 48271,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 73306, 0, 3,
                                                                       71206, 46921, 71416,
                                                                       26296, 26446, 48421,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 73621, 0, 3,
                                                                       71416, 47071, 71626,
                                                                       26446, 26596, 48646,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 73936, 0, 3,
                                                                       71626, 47221, 71836,
                                                                       26596, 26746, 48871,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 74251, 0, 3,
                                                                       71836, 47371, 72046,
                                                                       26746, 26896, 49096,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 74566, 0, 3,
                                                                       72046, 47521, 72256,
                                                                       26896, 27046, 49321,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 74881, 0, 3,
                                                                       72256, 47671, 72466,
                                                                       27046, 27196, 49546,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 75196, 0, 3,
                                                                       72466, 47821, 72676,
                                                                       27196, 27346, 49771,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 75511, 0, 3,
                                                                       72676, 47971, 72886,
                                                                       27346, 27496, 49996,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 75826, 0, 3,
                                                                       72886, 48121, 73096,
                                                                       27496, 27646, 50221,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 76141, 0, 3,
                                                                       73306, 48421, 73621,
                                                                       27946, 28156, 50446,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 76582, 0, 3,
                                                                       73621, 48646, 73936,
                                                                       28156, 28366, 50761,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 77023, 0, 3,
                                                                       73936, 48871, 74251,
                                                                       28366, 28576, 51076,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 77464, 0, 3,
                                                                       74251, 49096, 74566,
                                                                       28576, 28786, 51391,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 77905, 0, 3,
                                                                       74566, 49321, 74881,
                                                                       28786, 28996, 51706,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 78346, 0, 3,
                                                                       74881, 49546, 75196,
                                                                       28996, 29206, 52021,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 78787, 0, 3,
                                                                       75196, 49771, 75511,
                                                                       29206, 29416, 52336,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 79228, 0, 3,
                                                                       75511, 49996, 75826,
                                                                       29416, 29626, 52651,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 79669, 0, 3,
                                                                       76141, 50446, 76582,
                                                                       30046, 30326, 52966,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 80257, 0, 3,
                                                                       76582, 50761, 77023,
                                                                       30326, 30606, 53386,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 80845, 0, 3,
                                                                       77023, 51076, 77464,
                                                                       30606, 30886, 53806,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 81433, 0, 3,
                                                                       77464, 51391, 77905,
                                                                       30886, 31166, 54226,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 82021, 0, 3,
                                                                       77905, 51706, 78346,
                                                                       31166, 31446, 54646,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 82609, 0, 3,
                                                                       78346, 52021, 78787,
                                                                       31446, 31726, 55066,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 83197, 0, 3,
                                                                       78787, 52336, 79228,
                                                                       31726, 32006, 55486,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 83785, 0, 3,
                                                                       79669, 52966, 80257,
                                                                       32566, 32926, 55906,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 84541, 0, 3,
                                                                       80257, 53386, 80845,
                                                                       32926, 33286, 56446,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 85297, 0, 3,
                                                                       80845, 53806, 81433,
                                                                       33286, 33646, 56986,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 86053, 0, 3,
                                                                       81433, 54226, 82021,
                                                                       33646, 34006, 57526,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 86809, 0, 3,
                                                                       82021, 54646, 82609,
                                                                       34006, 34366, 58066,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 87565, 0, 3,
                                                                       82609, 55066, 83197,
                                                                       34366, 34726, 58606,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 88321, 0, 3,
                                                                       83785, 55906, 84541,
                                                                       35446, 35896, 59146,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 89266, 0, 3,
                                                                       84541, 56446, 85297,
                                                                       35896, 36346, 59821,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 90211, 0, 3,
                                                                       85297, 56986, 86053,
                                                                       36346, 36796, 60496,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 91156, 0, 3,
                                                                       86053, 57526, 86809,
                                                                       36796, 37246, 61171,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 92101, 0, 3,
                                                                       86809, 58066, 87565,
                                                                       37246, 37696, 61846,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 93046, 0, 3,
                                                                       88321, 59146, 89266,
                                                                       38596, 39146, 62521,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 94201, 0, 3,
                                                                       89266, 59821, 90211,
                                                                       39146, 39696, 63346,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 95356, 0, 3,
                                                                       90211, 60496, 91156,
                                                                       39696, 40246, 64171,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 96511, 0, 3,
                                                                       91156, 61171, 92101,
                                                                       40246, 40796, 64996,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 97666, 0, 3,
                                                                       93046, 62521, 94201,
                                                                       41896, 42556, 65821,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 99052, 0, 3,
                                                                       94201, 63346, 95356,
                                                                       42556, 43216, 66811,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 100438, 0, 3,
                                                                       95356, 64171, 96511,
                                                                       43216, 43876, 67801,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101824, 3, 45196,
                                                                       45211, 68833, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101852, 3, 45211,
                                                                       45226, 68854, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101880, 3, 45226,
                                                                       45241, 68875, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101908, 3, 45241,
                                                                       45256, 68896, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101936, 3, 45256,
                                                                       45271, 68917, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101964, 3, 45271,
                                                                       45286, 68938, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101992, 3, 45286,
                                                                       45301, 68959, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102020, 3, 45301,
                                                                       45316, 68980, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102048, 3, 45316,
                                                                       45331, 69001, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102076, 3, 45331,
                                                                       45346, 69022, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102104, 3, 45346,
                                                                       45361, 69043, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102132, 0, 3,
                                                                       101824, 68833, 101852,
                                                                       45391, 45436, 69190,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102216, 0, 3,
                                                                       101852, 68854, 101880,
                                                                       45436, 45481, 69253,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102300, 0, 3,
                                                                       101880, 68875, 101908,
                                                                       45481, 45526, 69316,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102384, 0, 3,
                                                                       101908, 68896, 101936,
                                                                       45526, 45571, 69379,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102468, 0, 3,
                                                                       101936, 68917, 101964,
                                                                       45571, 45616, 69442,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102552, 0, 3,
                                                                       101964, 68938, 101992,
                                                                       45616, 45661, 69505,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102636, 0, 3,
                                                                       101992, 68959, 102020,
                                                                       45661, 45706, 69568,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102720, 0, 3,
                                                                       102020, 68980, 102048,
                                                                       45706, 45751, 69631,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102804, 0, 3,
                                                                       102048, 69001, 102076,
                                                                       45751, 45796, 69694,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102888, 0, 3,
                                                                       102076, 69022, 102104,
                                                                       45796, 45841, 69757,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 102972, 0, 3,
                                                                       102132, 69190, 102216,
                                                                       45931, 46021, 70072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103140, 0, 3,
                                                                       102216, 69253, 102300,
                                                                       46021, 46111, 70198,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103308, 0, 3,
                                                                       102300, 69316, 102384,
                                                                       46111, 46201, 70324,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103476, 0, 3,
                                                                       102384, 69379, 102468,
                                                                       46201, 46291, 70450,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103644, 0, 3,
                                                                       102468, 69442, 102552,
                                                                       46291, 46381, 70576,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103812, 0, 3,
                                                                       102552, 69505, 102636,
                                                                       46381, 46471, 70702,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103980, 0, 3,
                                                                       102636, 69568, 102720,
                                                                       46471, 46561, 70828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 104148, 0, 3,
                                                                       102720, 69631, 102804,
                                                                       46561, 46651, 70954,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 104316, 0, 3,
                                                                       102804, 69694, 102888,
                                                                       46651, 46741, 71080,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 104484, 0, 3,
                                                                       102972, 70072, 103140,
                                                                       46921, 47071, 71626,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 104764, 0, 3,
                                                                       103140, 70198, 103308,
                                                                       47071, 47221, 71836,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 105044, 0, 3,
                                                                       103308, 70324, 103476,
                                                                       47221, 47371, 72046,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 105324, 0, 3,
                                                                       103476, 70450, 103644,
                                                                       47371, 47521, 72256,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 105604, 0, 3,
                                                                       103644, 70576, 103812,
                                                                       47521, 47671, 72466,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 105884, 0, 3,
                                                                       103812, 70702, 103980,
                                                                       47671, 47821, 72676,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 106164, 0, 3,
                                                                       103980, 70828, 104148,
                                                                       47821, 47971, 72886,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 106444, 0, 3,
                                                                       104148, 70954, 104316,
                                                                       47971, 48121, 73096,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 106724, 0, 3,
                                                                       104484, 71626, 104764,
                                                                       48421, 48646, 73936,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 107144, 0, 3,
                                                                       104764, 71836, 105044,
                                                                       48646, 48871, 74251,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 107564, 0, 3,
                                                                       105044, 72046, 105324,
                                                                       48871, 49096, 74566,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 107984, 0, 3,
                                                                       105324, 72256, 105604,
                                                                       49096, 49321, 74881,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 108404, 0, 3,
                                                                       105604, 72466, 105884,
                                                                       49321, 49546, 75196,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 108824, 0, 3,
                                                                       105884, 72676, 106164,
                                                                       49546, 49771, 75511,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 109244, 0, 3,
                                                                       106164, 72886, 106444,
                                                                       49771, 49996, 75826,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 109664, 0, 3,
                                                                       106724, 73936, 107144,
                                                                       50446, 50761, 77023,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 110252, 0, 3,
                                                                       107144, 74251, 107564,
                                                                       50761, 51076, 77464,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 110840, 0, 3,
                                                                       107564, 74566, 107984,
                                                                       51076, 51391, 77905,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 111428, 0, 3,
                                                                       107984, 74881, 108404,
                                                                       51391, 51706, 78346,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 112016, 0, 3,
                                                                       108404, 75196, 108824,
                                                                       51706, 52021, 78787,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 112604, 0, 3,
                                                                       108824, 75511, 109244,
                                                                       52021, 52336, 79228,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 113192, 0, 3,
                                                                       109664, 77023, 110252,
                                                                       52966, 53386, 80845,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 113976, 0, 3,
                                                                       110252, 77464, 110840,
                                                                       53386, 53806, 81433,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 114760, 0, 3,
                                                                       110840, 77905, 111428,
                                                                       53806, 54226, 82021,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 115544, 0, 3,
                                                                       111428, 78346, 112016,
                                                                       54226, 54646, 82609,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 116328, 0, 3,
                                                                       112016, 78787, 112604,
                                                                       54646, 55066, 83197,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 117112, 0, 3,
                                                                       113192, 80845, 113976,
                                                                       55906, 56446, 85297,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 118120, 0, 3,
                                                                       113976, 81433, 114760,
                                                                       56446, 56986, 86053,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 119128, 0, 3,
                                                                       114760, 82021, 115544,
                                                                       56986, 57526, 86809,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 120136, 0, 3,
                                                                       115544, 82609, 116328,
                                                                       57526, 58066, 87565,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 121144, 0, 3,
                                                                       117112, 85297, 118120,
                                                                       59146, 59821, 90211,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 122404, 0, 3,
                                                                       118120, 86053, 119128,
                                                                       59821, 60496, 91156,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 123664, 0, 3,
                                                                       119128, 86809, 120136,
                                                                       60496, 61171, 92101,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 124924, 0, 3,
                                                                       121144, 90211, 122404,
                                                                       62521, 63346, 95356,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 126464, 0, 3,
                                                                       122404, 91156, 123664,
                                                                       63346, 64171, 96511,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 128004, 0, 3,
                                                                       124924, 95356, 126464,
                                                                       65821, 66811, 100438,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129852, 3, 68791,
                                                                       68812, 101824, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129888, 3, 68812,
                                                                       68833, 101852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129924, 3, 68833,
                                                                       68854, 101880, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129960, 3, 68854,
                                                                       68875, 101908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129996, 3, 68875,
                                                                       68896, 101936, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130032, 3, 68896,
                                                                       68917, 101964, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130068, 3, 68917,
                                                                       68938, 101992, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130104, 3, 68938,
                                                                       68959, 102020, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130140, 3, 68959,
                                                                       68980, 102048, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130176, 3, 68980,
                                                                       69001, 102076, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130212, 3, 69001,
                                                                       69022, 102104, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130248, 0, 3,
                                                                       129852, 101824, 129888,
                                                                       69064, 69127, 102132,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130356, 0, 3,
                                                                       129888, 101852, 129924,
                                                                       69127, 69190, 102216,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130464, 0, 3,
                                                                       129924, 101880, 129960,
                                                                       69190, 69253, 102300,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130572, 0, 3,
                                                                       129960, 101908, 129996,
                                                                       69253, 69316, 102384,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130680, 0, 3,
                                                                       129996, 101936, 130032,
                                                                       69316, 69379, 102468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130788, 0, 3,
                                                                       130032, 101964, 130068,
                                                                       69379, 69442, 102552,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130896, 0, 3,
                                                                       130068, 101992, 130104,
                                                                       69442, 69505, 102636,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 131004, 0, 3,
                                                                       130104, 102020, 130140,
                                                                       69505, 69568, 102720,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 131112, 0, 3,
                                                                       130140, 102048, 130176,
                                                                       69568, 69631, 102804,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 131220, 0, 3,
                                                                       130176, 102076, 130212,
                                                                       69631, 69694, 102888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 131328, 0, 3,
                                                                       130248, 102132, 130356,
                                                                       69820, 69946, 102972,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 131544, 0, 3,
                                                                       130356, 102216, 130464,
                                                                       69946, 70072, 103140,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 131760, 0, 3,
                                                                       130464, 102300, 130572,
                                                                       70072, 70198, 103308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 131976, 0, 3,
                                                                       130572, 102384, 130680,
                                                                       70198, 70324, 103476,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 132192, 0, 3,
                                                                       130680, 102468, 130788,
                                                                       70324, 70450, 103644,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 132408, 0, 3,
                                                                       130788, 102552, 130896,
                                                                       70450, 70576, 103812,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 132624, 0, 3,
                                                                       130896, 102636, 131004,
                                                                       70576, 70702, 103980,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 132840, 0, 3,
                                                                       131004, 102720, 131112,
                                                                       70702, 70828, 104148,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 133056, 0, 3,
                                                                       131112, 102804, 131220,
                                                                       70828, 70954, 104316,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 133272, 0, 3,
                                                                       131328, 102972, 131544,
                                                                       71206, 71416, 104484,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 133632, 0, 3,
                                                                       131544, 103140, 131760,
                                                                       71416, 71626, 104764,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 133992, 0, 3,
                                                                       131760, 103308, 131976,
                                                                       71626, 71836, 105044,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 134352, 0, 3,
                                                                       131976, 103476, 132192,
                                                                       71836, 72046, 105324,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 134712, 0, 3,
                                                                       132192, 103644, 132408,
                                                                       72046, 72256, 105604,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 135072, 0, 3,
                                                                       132408, 103812, 132624,
                                                                       72256, 72466, 105884,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 135432, 0, 3,
                                                                       132624, 103980, 132840,
                                                                       72466, 72676, 106164,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 135792, 0, 3,
                                                                       132840, 104148, 133056,
                                                                       72676, 72886, 106444,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 136152, 0, 3,
                                                                       133272, 104484, 133632,
                                                                       73306, 73621, 106724,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 136692, 0, 3,
                                                                       133632, 104764, 133992,
                                                                       73621, 73936, 107144,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 137232, 0, 3,
                                                                       133992, 105044, 134352,
                                                                       73936, 74251, 107564,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 137772, 0, 3,
                                                                       134352, 105324, 134712,
                                                                       74251, 74566, 107984,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 138312, 0, 3,
                                                                       134712, 105604, 135072,
                                                                       74566, 74881, 108404,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 138852, 0, 3,
                                                                       135072, 105884, 135432,
                                                                       74881, 75196, 108824,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 139392, 0, 3,
                                                                       135432, 106164, 135792,
                                                                       75196, 75511, 109244,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 139932, 0, 3,
                                                                       136152, 106724, 136692,
                                                                       76141, 76582, 109664,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 140688, 0, 3,
                                                                       136692, 107144, 137232,
                                                                       76582, 77023, 110252,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 141444, 0, 3,
                                                                       137232, 107564, 137772,
                                                                       77023, 77464, 110840,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 142200, 0, 3,
                                                                       137772, 107984, 138312,
                                                                       77464, 77905, 111428,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 142956, 0, 3,
                                                                       138312, 108404, 138852,
                                                                       77905, 78346, 112016,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 143712, 0, 3,
                                                                       138852, 108824, 139392,
                                                                       78346, 78787, 112604,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 144468, 0, 3,
                                                                       139932, 109664, 140688,
                                                                       79669, 80257, 113192,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 145476, 0, 3,
                                                                       140688, 110252, 141444,
                                                                       80257, 80845, 113976,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 146484, 0, 3,
                                                                       141444, 110840, 142200,
                                                                       80845, 81433, 114760,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 147492, 0, 3,
                                                                       142200, 111428, 142956,
                                                                       81433, 82021, 115544,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 148500, 0, 3,
                                                                       142956, 112016, 143712,
                                                                       82021, 82609, 116328,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 149508, 0, 3,
                                                                       144468, 113192, 145476,
                                                                       83785, 84541, 117112,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 150804, 0, 3,
                                                                       145476, 113976, 146484,
                                                                       84541, 85297, 118120,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 152100, 0, 3,
                                                                       146484, 114760, 147492,
                                                                       85297, 86053, 119128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 153396, 0, 3,
                                                                       147492, 115544, 148500,
                                                                       86053, 86809, 120136,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 154692, 0, 3,
                                                                       149508, 117112, 150804,
                                                                       88321, 89266, 121144,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 156312, 0, 3,
                                                                       150804, 118120, 152100,
                                                                       89266, 90211, 122404,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 157932, 0, 3,
                                                                       152100, 119128, 153396,
                                                                       90211, 91156, 123664,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 159552, 0, 3,
                                                                       154692, 121144, 156312,
                                                                       93046, 94201, 124924,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 161532, 0, 3,
                                                                       156312, 122404, 157932,
                                                                       94201, 95356, 126464,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 163512, 0, 3,
                                                                       159552, 124924, 161532,
                                                                       97666, 99052, 128004,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 165888, 139932, 756, ncols);

                    simdfunc::contract_primitives(buffer, 166959, 144468, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 168387, 149508, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 170223, 154692, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 172518, 159552, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 175323, 163512, 2376, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 166644, 165888, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 167967, 166959, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 169683, 168387, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 171843, 170223, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 174498, 172518, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 177699, 175323, 66, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 178689, 166644, 167967, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 179634, 167967, 169683, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 180894, 169683, 171843, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 182514, 171843, 174498, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 184539, 174498, 177699, 15,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 187014, 178689, 179634, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 188904, 179634, 180894, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 191424, 180894, 182514, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 194664, 182514, 184539, 15,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 198714, 187014, 188904, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 201864, 188904, 191424, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 206064, 191424, 194664, 15,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 211464, 198714, 201864, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 216189, 201864, 206064, 15,
                                             nmax);

        simdtrf::compute_hrr_hh_out_of_first(buffer, coordinates, 222489, 211464, 216189, 15,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 229104, 222489, 21, 15, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 229104, 165, nmax);
    }

    for (size_t m = 0; m < 1815; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
