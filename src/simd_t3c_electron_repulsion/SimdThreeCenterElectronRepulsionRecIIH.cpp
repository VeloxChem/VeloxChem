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


#include "SimdThreeCenterElectronRepulsionRecIIH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
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
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iih_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iih_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 225317, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 225317, 123768, 11767, dimensions);

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

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2725, 0, 3, 1823,
                                                                       1878, 2263, 2329, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2803, 0, 3, 1878,
                                                                       1933, 2329, 2395, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2881, 0, 3, 1933,
                                                                       1988, 2395, 2461, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2959, 0, 3, 1988,
                                                                       2043, 2461, 2527, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3037, 0, 3, 2043,
                                                                       2098, 2527, 2593, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3115, 0, 3, 2098,
                                                                       2153, 2593, 2659, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3193, 0, 3, 2263,
                                                                       2329, 2725, 2803, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3284, 0, 3, 2329,
                                                                       2395, 2803, 2881, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3375, 0, 3, 2395,
                                                                       2461, 2881, 2959, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3466, 0, 3, 2461,
                                                                       2527, 2959, 3037, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3557, 0, 3, 2527,
                                                                       2593, 3037, 3115, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3648, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3651, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3654, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3657, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3660, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3663, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3666, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3669, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3672, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3675, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3678, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3681, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3684, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3687, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3690, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3693, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3696, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3699, 3, 8, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3708, 3, 9, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3717, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3726, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3735, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3744, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3753, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3762, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3771, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3780, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3789, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3798, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3807, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3816, 3, 21, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3825, 3, 22, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3834, 3, 23, 70,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3843, 3, 25, 73,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3861, 3, 28, 79,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3879, 3, 31, 85,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3897, 3, 34, 91,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3915, 3, 37, 97,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3933, 3, 40, 103,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3951, 3, 43, 109,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3969, 3, 46, 115,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3987, 3, 49, 121,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4005, 3, 52, 127,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4023, 3, 55, 133,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4041, 3, 58, 139,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4059, 3, 61, 145,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4077, 3, 64, 151,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4095, 3, 67, 157,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4113, 3, 73, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4143, 3, 79, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4173, 3, 85, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4203, 3, 91, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4233, 3, 97, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4263, 3, 103, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4293, 3, 109, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4323, 3, 115, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4353, 3, 121, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4383, 3, 127, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4413, 3, 133, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4443, 3, 139, 273,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4473, 3, 145, 283,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4503, 3, 151, 293,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4533, 3, 163, 303,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4578, 3, 173, 318,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4623, 3, 183, 333,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4668, 3, 193, 348,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4713, 3, 203, 363,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4758, 3, 213, 378,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4803, 3, 223, 393,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4848, 3, 233, 408,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4893, 3, 243, 423,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4938, 3, 253, 438,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4983, 3, 263, 453,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5028, 3, 273, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5073, 3, 283, 483,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5118, 3, 303, 498,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5181, 3, 318, 519,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5244, 3, 333, 540,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5307, 3, 348, 561,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5370, 3, 363, 582,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5433, 3, 378, 603,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5496, 3, 393, 624,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5559, 3, 408, 645,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5622, 3, 423, 666,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5685, 3, 438, 687,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5748, 3, 453, 708,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5811, 3, 468, 729,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5874, 3, 498, 750,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5958, 3, 519, 778,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6042, 3, 540, 806,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6126, 3, 561, 834,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6210, 3, 582, 862,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6294, 3, 603, 890,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6378, 3, 624, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6462, 3, 645, 946,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6546, 3, 666, 974,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6630, 3, 687,
                                                                       1002, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6714, 3, 708,
                                                                       1030, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6798, 3, 750,
                                                                       1058, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6906, 3, 778,
                                                                       1094, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7014, 3, 806,
                                                                       1130, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7122, 3, 834,
                                                                       1166, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7230, 3, 862,
                                                                       1202, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7338, 3, 890,
                                                                       1238, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7446, 3, 918,
                                                                       1274, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7554, 3, 946,
                                                                       1310, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7662, 3, 974,
                                                                       1346, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7770, 3, 1002,
                                                                       1382, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7878, 3, 1058,
                                                                       1418, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8013, 3, 1094,
                                                                       1463, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8148, 3, 1130,
                                                                       1508, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8283, 3, 1166,
                                                                       1553, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8418, 3, 1202,
                                                                       1598, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8553, 3, 1238,
                                                                       1643, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8688, 3, 1274,
                                                                       1688, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8823, 3, 1310,
                                                                       1733, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8958, 3, 1346,
                                                                       1778, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9093, 3, 1418,
                                                                       1823, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9258, 3, 1463,
                                                                       1878, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9423, 3, 1508,
                                                                       1933, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9588, 3, 1553,
                                                                       1988, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9753, 3, 1598,
                                                                       2043, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9918, 3, 1643,
                                                                       2098, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10083, 3, 1688,
                                                                       2153, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10248, 3, 1733,
                                                                       2208, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10413, 3, 1823,
                                                                       2263, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10611, 3, 1878,
                                                                       2329, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10809, 3, 1933,
                                                                       2395, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11007, 3, 1988,
                                                                       2461, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11205, 3, 2043,
                                                                       2527, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11403, 3, 2098,
                                                                       2593, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11601, 3, 2153,
                                                                       2659, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 11799, 3, 2263,
                                                                       2725, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12033, 3, 2329,
                                                                       2803, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12267, 3, 2395,
                                                                       2881, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12501, 3, 2461,
                                                                       2959, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12735, 3, 2527,
                                                                       3037, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12969, 3, 2593,
                                                                       3115, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 13203, 3, 2725,
                                                                       3193, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 13476, 3, 2803,
                                                                       3284, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 13749, 3, 2881,
                                                                       3375, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14022, 3, 2959,
                                                                       3466, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14295, 3, 3037,
                                                                       3557, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14568, 3, 8, 9,
                                                                       3654, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14574, 3, 9, 10,
                                                                       3657, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14580, 3, 10, 11,
                                                                       3660, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14586, 3, 11, 12,
                                                                       3663, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14592, 3, 12, 13,
                                                                       3666, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14598, 3, 13, 14,
                                                                       3669, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14604, 3, 14, 15,
                                                                       3672, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14610, 3, 15, 16,
                                                                       3675, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14616, 3, 16, 17,
                                                                       3678, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14622, 3, 17, 18,
                                                                       3681, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14628, 3, 18, 19,
                                                                       3684, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14634, 3, 19, 20,
                                                                       3687, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14640, 3, 20, 21,
                                                                       3690, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14646, 3, 21, 22,
                                                                       3693, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14652, 3, 22, 23,
                                                                       3696, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14658, 0, 3,
                                                                       14568, 3654, 14574, 3717,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14676, 0, 3,
                                                                       14574, 3657, 14580, 3726,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14694, 0, 3,
                                                                       14580, 3660, 14586, 3735,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14712, 0, 3,
                                                                       14586, 3663, 14592, 3744,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14730, 0, 3,
                                                                       14592, 3666, 14598, 3753,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14748, 0, 3,
                                                                       14598, 3669, 14604, 3762,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14766, 0, 3,
                                                                       14604, 3672, 14610, 3771,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14784, 0, 3,
                                                                       14610, 3675, 14616, 3780,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14802, 0, 3,
                                                                       14616, 3678, 14622, 3789,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14820, 0, 3,
                                                                       14622, 3681, 14628, 3798,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14838, 0, 3,
                                                                       14628, 3684, 14634, 3807,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14856, 0, 3,
                                                                       14634, 3687, 14640, 3816,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14874, 0, 3,
                                                                       14640, 3690, 14646, 3825,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14892, 0, 3,
                                                                       14646, 3693, 14652, 3834,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14910, 0, 3,
                                                                       14658, 3717, 14676, 73,
                                                                       79, 3879, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14946, 0, 3,
                                                                       14676, 3726, 14694, 79,
                                                                       85, 3897, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14982, 0, 3,
                                                                       14694, 3735, 14712, 85,
                                                                       91, 3915, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15018, 0, 3,
                                                                       14712, 3744, 14730, 91,
                                                                       97, 3933, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15054, 0, 3,
                                                                       14730, 3753, 14748, 97,
                                                                       103, 3951, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15090, 0, 3,
                                                                       14748, 3762, 14766, 103,
                                                                       109, 3969, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15126, 0, 3,
                                                                       14766, 3771, 14784, 109,
                                                                       115, 3987, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15162, 0, 3,
                                                                       14784, 3780, 14802, 115,
                                                                       121, 4005, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15198, 0, 3,
                                                                       14802, 3789, 14820, 121,
                                                                       127, 4023, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15234, 0, 3,
                                                                       14820, 3798, 14838, 127,
                                                                       133, 4041, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15270, 0, 3,
                                                                       14838, 3807, 14856, 133,
                                                                       139, 4059, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15306, 0, 3,
                                                                       14856, 3816, 14874, 139,
                                                                       145, 4077, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15342, 0, 3,
                                                                       14874, 3825, 14892, 145,
                                                                       151, 4095, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15378, 0, 3,
                                                                       14910, 3879, 14946, 163,
                                                                       173, 4173, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15438, 0, 3,
                                                                       14946, 3897, 14982, 173,
                                                                       183, 4203, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15498, 0, 3,
                                                                       14982, 3915, 15018, 183,
                                                                       193, 4233, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15558, 0, 3,
                                                                       15018, 3933, 15054, 193,
                                                                       203, 4263, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15618, 0, 3,
                                                                       15054, 3951, 15090, 203,
                                                                       213, 4293, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15678, 0, 3,
                                                                       15090, 3969, 15126, 213,
                                                                       223, 4323, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15738, 0, 3,
                                                                       15126, 3987, 15162, 223,
                                                                       233, 4353, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15798, 0, 3,
                                                                       15162, 4005, 15198, 233,
                                                                       243, 4383, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15858, 0, 3,
                                                                       15198, 4023, 15234, 243,
                                                                       253, 4413, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15918, 0, 3,
                                                                       15234, 4041, 15270, 253,
                                                                       263, 4443, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15978, 0, 3,
                                                                       15270, 4059, 15306, 263,
                                                                       273, 4473, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16038, 0, 3,
                                                                       15306, 4077, 15342, 273,
                                                                       283, 4503, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16098, 0, 3,
                                                                       15378, 4173, 15438, 303,
                                                                       318, 4623, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16188, 0, 3,
                                                                       15438, 4203, 15498, 318,
                                                                       333, 4668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16278, 0, 3,
                                                                       15498, 4233, 15558, 333,
                                                                       348, 4713, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16368, 0, 3,
                                                                       15558, 4263, 15618, 348,
                                                                       363, 4758, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16458, 0, 3,
                                                                       15618, 4293, 15678, 363,
                                                                       378, 4803, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16548, 0, 3,
                                                                       15678, 4323, 15738, 378,
                                                                       393, 4848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16638, 0, 3,
                                                                       15738, 4353, 15798, 393,
                                                                       408, 4893, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16728, 0, 3,
                                                                       15798, 4383, 15858, 408,
                                                                       423, 4938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16818, 0, 3,
                                                                       15858, 4413, 15918, 423,
                                                                       438, 4983, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16908, 0, 3,
                                                                       15918, 4443, 15978, 438,
                                                                       453, 5028, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16998, 0, 3,
                                                                       15978, 4473, 16038, 453,
                                                                       468, 5073, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17088, 0, 3,
                                                                       16098, 4623, 16188, 498,
                                                                       519, 5244, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17214, 0, 3,
                                                                       16188, 4668, 16278, 519,
                                                                       540, 5307, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17340, 0, 3,
                                                                       16278, 4713, 16368, 540,
                                                                       561, 5370, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17466, 0, 3,
                                                                       16368, 4758, 16458, 561,
                                                                       582, 5433, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17592, 0, 3,
                                                                       16458, 4803, 16548, 582,
                                                                       603, 5496, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17718, 0, 3,
                                                                       16548, 4848, 16638, 603,
                                                                       624, 5559, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17844, 0, 3,
                                                                       16638, 4893, 16728, 624,
                                                                       645, 5622, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17970, 0, 3,
                                                                       16728, 4938, 16818, 645,
                                                                       666, 5685, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18096, 0, 3,
                                                                       16818, 4983, 16908, 666,
                                                                       687, 5748, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18222, 0, 3,
                                                                       16908, 5028, 16998, 687,
                                                                       708, 5811, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18348, 0, 3,
                                                                       17088, 5244, 17214, 750,
                                                                       778, 6042, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18516, 0, 3,
                                                                       17214, 5307, 17340, 778,
                                                                       806, 6126, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18684, 0, 3,
                                                                       17340, 5370, 17466, 806,
                                                                       834, 6210, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18852, 0, 3,
                                                                       17466, 5433, 17592, 834,
                                                                       862, 6294, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19020, 0, 3,
                                                                       17592, 5496, 17718, 862,
                                                                       890, 6378, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19188, 0, 3,
                                                                       17718, 5559, 17844, 890,
                                                                       918, 6462, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19356, 0, 3,
                                                                       17844, 5622, 17970, 918,
                                                                       946, 6546, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19524, 0, 3,
                                                                       17970, 5685, 18096, 946,
                                                                       974, 6630, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19692, 0, 3,
                                                                       18096, 5748, 18222, 974,
                                                                       1002, 6714, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19860, 0, 3,
                                                                       18348, 6042, 18516, 1058,
                                                                       1094, 7014, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20076, 0, 3,
                                                                       18516, 6126, 18684, 1094,
                                                                       1130, 7122, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20292, 0, 3,
                                                                       18684, 6210, 18852, 1130,
                                                                       1166, 7230, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20508, 0, 3,
                                                                       18852, 6294, 19020, 1166,
                                                                       1202, 7338, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20724, 0, 3,
                                                                       19020, 6378, 19188, 1202,
                                                                       1238, 7446, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20940, 0, 3,
                                                                       19188, 6462, 19356, 1238,
                                                                       1274, 7554, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21156, 0, 3,
                                                                       19356, 6546, 19524, 1274,
                                                                       1310, 7662, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21372, 0, 3,
                                                                       19524, 6630, 19692, 1310,
                                                                       1346, 7770, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21588, 0, 3,
                                                                       19860, 7014, 20076, 1418,
                                                                       1463, 8148, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21858, 0, 3,
                                                                       20076, 7122, 20292, 1463,
                                                                       1508, 8283, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22128, 0, 3,
                                                                       20292, 7230, 20508, 1508,
                                                                       1553, 8418, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22398, 0, 3,
                                                                       20508, 7338, 20724, 1553,
                                                                       1598, 8553, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22668, 0, 3,
                                                                       20724, 7446, 20940, 1598,
                                                                       1643, 8688, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22938, 0, 3,
                                                                       20940, 7554, 21156, 1643,
                                                                       1688, 8823, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23208, 0, 3,
                                                                       21156, 7662, 21372, 1688,
                                                                       1733, 8958, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23478, 0, 3,
                                                                       21588, 8148, 21858, 1823,
                                                                       1878, 9423, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23808, 0, 3,
                                                                       21858, 8283, 22128, 1878,
                                                                       1933, 9588, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24138, 0, 3,
                                                                       22128, 8418, 22398, 1933,
                                                                       1988, 9753, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24468, 0, 3,
                                                                       22398, 8553, 22668, 1988,
                                                                       2043, 9918, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24798, 0, 3,
                                                                       22668, 8688, 22938, 2043,
                                                                       2098, 10083, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25128, 0, 3,
                                                                       22938, 8823, 23208, 2098,
                                                                       2153, 10248, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 25458, 0, 3,
                                                                       23478, 9423, 23808, 2263,
                                                                       2329, 10809, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 25854, 0, 3,
                                                                       23808, 9588, 24138, 2329,
                                                                       2395, 11007, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26250, 0, 3,
                                                                       24138, 9753, 24468, 2395,
                                                                       2461, 11205, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26646, 0, 3,
                                                                       24468, 9918, 24798, 2461,
                                                                       2527, 11403, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27042, 0, 3,
                                                                       24798, 10083, 25128, 2527,
                                                                       2593, 11601, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 27438, 0, 3,
                                                                       25458, 10809, 25854, 2725,
                                                                       2803, 12267, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 27906, 0, 3,
                                                                       25854, 11007, 26250, 2803,
                                                                       2881, 12501, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 28374, 0, 3,
                                                                       26250, 11205, 26646, 2881,
                                                                       2959, 12735, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 28842, 0, 3,
                                                                       26646, 11403, 27042, 2959,
                                                                       3037, 12969, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 29310, 0, 3,
                                                                       27438, 12267, 27906, 3193,
                                                                       3284, 13749, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 29856, 0, 3,
                                                                       27906, 12501, 28374, 3284,
                                                                       3375, 14022, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 30402, 0, 3,
                                                                       28374, 12735, 28842, 3375,
                                                                       3466, 14295, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30948, 3, 3648,
                                                                       3651, 14568, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30958, 3, 3651,
                                                                       3654, 14574, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30968, 3, 3654,
                                                                       3657, 14580, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30978, 3, 3657,
                                                                       3660, 14586, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30988, 3, 3660,
                                                                       3663, 14592, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30998, 3, 3663,
                                                                       3666, 14598, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31008, 3, 3666,
                                                                       3669, 14604, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31018, 3, 3669,
                                                                       3672, 14610, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31028, 3, 3672,
                                                                       3675, 14616, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31038, 3, 3675,
                                                                       3678, 14622, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31048, 3, 3678,
                                                                       3681, 14628, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31058, 3, 3681,
                                                                       3684, 14634, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31068, 3, 3684,
                                                                       3687, 14640, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31078, 3, 3687,
                                                                       3690, 14646, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31088, 3, 3690,
                                                                       3693, 14652, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31098, 0, 3,
                                                                       30948, 14568, 30958, 3699,
                                                                       3708, 14658, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31128, 0, 3,
                                                                       30958, 14574, 30968, 3708,
                                                                       3717, 14676, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31158, 0, 3,
                                                                       30968, 14580, 30978, 3717,
                                                                       3726, 14694, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31188, 0, 3,
                                                                       30978, 14586, 30988, 3726,
                                                                       3735, 14712, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31218, 0, 3,
                                                                       30988, 14592, 30998, 3735,
                                                                       3744, 14730, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31248, 0, 3,
                                                                       30998, 14598, 31008, 3744,
                                                                       3753, 14748, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31278, 0, 3,
                                                                       31008, 14604, 31018, 3753,
                                                                       3762, 14766, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31308, 0, 3,
                                                                       31018, 14610, 31028, 3762,
                                                                       3771, 14784, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31338, 0, 3,
                                                                       31028, 14616, 31038, 3771,
                                                                       3780, 14802, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31368, 0, 3,
                                                                       31038, 14622, 31048, 3780,
                                                                       3789, 14820, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31398, 0, 3,
                                                                       31048, 14628, 31058, 3789,
                                                                       3798, 14838, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31428, 0, 3,
                                                                       31058, 14634, 31068, 3798,
                                                                       3807, 14856, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31458, 0, 3,
                                                                       31068, 14640, 31078, 3807,
                                                                       3816, 14874, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31488, 0, 3,
                                                                       31078, 14646, 31088, 3816,
                                                                       3825, 14892, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31518, 0, 3,
                                                                       31098, 14658, 31128, 3843,
                                                                       3861, 14910, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31578, 0, 3,
                                                                       31128, 14676, 31158, 3861,
                                                                       3879, 14946, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31638, 0, 3,
                                                                       31158, 14694, 31188, 3879,
                                                                       3897, 14982, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31698, 0, 3,
                                                                       31188, 14712, 31218, 3897,
                                                                       3915, 15018, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31758, 0, 3,
                                                                       31218, 14730, 31248, 3915,
                                                                       3933, 15054, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31818, 0, 3,
                                                                       31248, 14748, 31278, 3933,
                                                                       3951, 15090, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31878, 0, 3,
                                                                       31278, 14766, 31308, 3951,
                                                                       3969, 15126, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31938, 0, 3,
                                                                       31308, 14784, 31338, 3969,
                                                                       3987, 15162, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31998, 0, 3,
                                                                       31338, 14802, 31368, 3987,
                                                                       4005, 15198, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32058, 0, 3,
                                                                       31368, 14820, 31398, 4005,
                                                                       4023, 15234, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32118, 0, 3,
                                                                       31398, 14838, 31428, 4023,
                                                                       4041, 15270, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32178, 0, 3,
                                                                       31428, 14856, 31458, 4041,
                                                                       4059, 15306, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32238, 0, 3,
                                                                       31458, 14874, 31488, 4059,
                                                                       4077, 15342, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32298, 0, 3,
                                                                       31518, 14910, 31578, 4113,
                                                                       4143, 15378, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32398, 0, 3,
                                                                       31578, 14946, 31638, 4143,
                                                                       4173, 15438, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32498, 0, 3,
                                                                       31638, 14982, 31698, 4173,
                                                                       4203, 15498, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32598, 0, 3,
                                                                       31698, 15018, 31758, 4203,
                                                                       4233, 15558, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32698, 0, 3,
                                                                       31758, 15054, 31818, 4233,
                                                                       4263, 15618, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32798, 0, 3,
                                                                       31818, 15090, 31878, 4263,
                                                                       4293, 15678, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32898, 0, 3,
                                                                       31878, 15126, 31938, 4293,
                                                                       4323, 15738, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32998, 0, 3,
                                                                       31938, 15162, 31998, 4323,
                                                                       4353, 15798, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33098, 0, 3,
                                                                       31998, 15198, 32058, 4353,
                                                                       4383, 15858, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33198, 0, 3,
                                                                       32058, 15234, 32118, 4383,
                                                                       4413, 15918, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33298, 0, 3,
                                                                       32118, 15270, 32178, 4413,
                                                                       4443, 15978, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33398, 0, 3,
                                                                       32178, 15306, 32238, 4443,
                                                                       4473, 16038, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 33498, 0, 3,
                                                                       32298, 15378, 32398, 4533,
                                                                       4578, 16098, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 33648, 0, 3,
                                                                       32398, 15438, 32498, 4578,
                                                                       4623, 16188, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 33798, 0, 3,
                                                                       32498, 15498, 32598, 4623,
                                                                       4668, 16278, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 33948, 0, 3,
                                                                       32598, 15558, 32698, 4668,
                                                                       4713, 16368, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34098, 0, 3,
                                                                       32698, 15618, 32798, 4713,
                                                                       4758, 16458, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34248, 0, 3,
                                                                       32798, 15678, 32898, 4758,
                                                                       4803, 16548, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34398, 0, 3,
                                                                       32898, 15738, 32998, 4803,
                                                                       4848, 16638, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34548, 0, 3,
                                                                       32998, 15798, 33098, 4848,
                                                                       4893, 16728, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34698, 0, 3,
                                                                       33098, 15858, 33198, 4893,
                                                                       4938, 16818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34848, 0, 3,
                                                                       33198, 15918, 33298, 4938,
                                                                       4983, 16908, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34998, 0, 3,
                                                                       33298, 15978, 33398, 4983,
                                                                       5028, 16998, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35148, 0, 3,
                                                                       33498, 16098, 33648, 5118,
                                                                       5181, 17088, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35358, 0, 3,
                                                                       33648, 16188, 33798, 5181,
                                                                       5244, 17214, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35568, 0, 3,
                                                                       33798, 16278, 33948, 5244,
                                                                       5307, 17340, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35778, 0, 3,
                                                                       33948, 16368, 34098, 5307,
                                                                       5370, 17466, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35988, 0, 3,
                                                                       34098, 16458, 34248, 5370,
                                                                       5433, 17592, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36198, 0, 3,
                                                                       34248, 16548, 34398, 5433,
                                                                       5496, 17718, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36408, 0, 3,
                                                                       34398, 16638, 34548, 5496,
                                                                       5559, 17844, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36618, 0, 3,
                                                                       34548, 16728, 34698, 5559,
                                                                       5622, 17970, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36828, 0, 3,
                                                                       34698, 16818, 34848, 5622,
                                                                       5685, 18096, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37038, 0, 3,
                                                                       34848, 16908, 34998, 5685,
                                                                       5748, 18222, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37248, 0, 3,
                                                                       35148, 17088, 35358, 5874,
                                                                       5958, 18348, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37528, 0, 3,
                                                                       35358, 17214, 35568, 5958,
                                                                       6042, 18516, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37808, 0, 3,
                                                                       35568, 17340, 35778, 6042,
                                                                       6126, 18684, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38088, 0, 3,
                                                                       35778, 17466, 35988, 6126,
                                                                       6210, 18852, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38368, 0, 3,
                                                                       35988, 17592, 36198, 6210,
                                                                       6294, 19020, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38648, 0, 3,
                                                                       36198, 17718, 36408, 6294,
                                                                       6378, 19188, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38928, 0, 3,
                                                                       36408, 17844, 36618, 6378,
                                                                       6462, 19356, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39208, 0, 3,
                                                                       36618, 17970, 36828, 6462,
                                                                       6546, 19524, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39488, 0, 3,
                                                                       36828, 18096, 37038, 6546,
                                                                       6630, 19692, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 39768, 0, 3,
                                                                       37248, 18348, 37528, 6798,
                                                                       6906, 19860, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40128, 0, 3,
                                                                       37528, 18516, 37808, 6906,
                                                                       7014, 20076, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40488, 0, 3,
                                                                       37808, 18684, 38088, 7014,
                                                                       7122, 20292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40848, 0, 3,
                                                                       38088, 18852, 38368, 7122,
                                                                       7230, 20508, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41208, 0, 3,
                                                                       38368, 19020, 38648, 7230,
                                                                       7338, 20724, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41568, 0, 3,
                                                                       38648, 19188, 38928, 7338,
                                                                       7446, 20940, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41928, 0, 3,
                                                                       38928, 19356, 39208, 7446,
                                                                       7554, 21156, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42288, 0, 3,
                                                                       39208, 19524, 39488, 7554,
                                                                       7662, 21372, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 42648, 0, 3,
                                                                       39768, 19860, 40128, 7878,
                                                                       8013, 21588, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 43098, 0, 3,
                                                                       40128, 20076, 40488, 8013,
                                                                       8148, 21858, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 43548, 0, 3,
                                                                       40488, 20292, 40848, 8148,
                                                                       8283, 22128, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 43998, 0, 3,
                                                                       40848, 20508, 41208, 8283,
                                                                       8418, 22398, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44448, 0, 3,
                                                                       41208, 20724, 41568, 8418,
                                                                       8553, 22668, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44898, 0, 3,
                                                                       41568, 20940, 41928, 8553,
                                                                       8688, 22938, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 45348, 0, 3,
                                                                       41928, 21156, 42288, 8688,
                                                                       8823, 23208, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 45798, 0, 3,
                                                                       42648, 21588, 43098, 9093,
                                                                       9258, 23478, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 46348, 0, 3,
                                                                       43098, 21858, 43548, 9258,
                                                                       9423, 23808, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 46898, 0, 3,
                                                                       43548, 22128, 43998, 9423,
                                                                       9588, 24138, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 47448, 0, 3,
                                                                       43998, 22398, 44448, 9588,
                                                                       9753, 24468, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 47998, 0, 3,
                                                                       44448, 22668, 44898, 9753,
                                                                       9918, 24798, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 48548, 0, 3,
                                                                       44898, 22938, 45348, 9918,
                                                                       10083, 25128, ncols,
                                                                       gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 49098, 0, 3,
                                                                       45798, 23478, 46348,
                                                                       10413, 10611, 25458,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 49758, 0, 3,
                                                                       46348, 23808, 46898,
                                                                       10611, 10809, 25854,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 50418, 0, 3,
                                                                       46898, 24138, 47448,
                                                                       10809, 11007, 26250,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 51078, 0, 3,
                                                                       47448, 24468, 47998,
                                                                       11007, 11205, 26646,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 51738, 0, 3,
                                                                       47998, 24798, 48548,
                                                                       11205, 11403, 27042,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 52398, 0, 3,
                                                                       49098, 25458, 49758,
                                                                       11799, 12033, 27438,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 53178, 0, 3,
                                                                       49758, 25854, 50418,
                                                                       12033, 12267, 27906,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 53958, 0, 3,
                                                                       50418, 26250, 51078,
                                                                       12267, 12501, 28374,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 54738, 0, 3,
                                                                       51078, 26646, 51738,
                                                                       12501, 12735, 28842,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 55518, 0, 3,
                                                                       52398, 27438, 53178,
                                                                       13203, 13476, 29310,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 56428, 0, 3,
                                                                       53178, 27906, 53958,
                                                                       13476, 13749, 29856,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 57338, 0, 3,
                                                                       53958, 28374, 54738,
                                                                       13749, 14022, 30402,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58248, 3, 14568,
                                                                       14574, 30968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58263, 3, 14574,
                                                                       14580, 30978, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58278, 3, 14580,
                                                                       14586, 30988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58293, 3, 14586,
                                                                       14592, 30998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58308, 3, 14592,
                                                                       14598, 31008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58323, 3, 14598,
                                                                       14604, 31018, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58338, 3, 14604,
                                                                       14610, 31028, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58353, 3, 14610,
                                                                       14616, 31038, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58368, 3, 14616,
                                                                       14622, 31048, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58383, 3, 14622,
                                                                       14628, 31058, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58398, 3, 14628,
                                                                       14634, 31068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58413, 3, 14634,
                                                                       14640, 31078, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58428, 3, 14640,
                                                                       14646, 31088, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58443, 0, 3,
                                                                       58248, 30968, 58263,
                                                                       14658, 14676, 31158,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58488, 0, 3,
                                                                       58263, 30978, 58278,
                                                                       14676, 14694, 31188,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58533, 0, 3,
                                                                       58278, 30988, 58293,
                                                                       14694, 14712, 31218,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58578, 0, 3,
                                                                       58293, 30998, 58308,
                                                                       14712, 14730, 31248,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58623, 0, 3,
                                                                       58308, 31008, 58323,
                                                                       14730, 14748, 31278,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58668, 0, 3,
                                                                       58323, 31018, 58338,
                                                                       14748, 14766, 31308,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58713, 0, 3,
                                                                       58338, 31028, 58353,
                                                                       14766, 14784, 31338,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58758, 0, 3,
                                                                       58353, 31038, 58368,
                                                                       14784, 14802, 31368,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58803, 0, 3,
                                                                       58368, 31048, 58383,
                                                                       14802, 14820, 31398,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58848, 0, 3,
                                                                       58383, 31058, 58398,
                                                                       14820, 14838, 31428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58893, 0, 3,
                                                                       58398, 31068, 58413,
                                                                       14838, 14856, 31458,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58938, 0, 3,
                                                                       58413, 31078, 58428,
                                                                       14856, 14874, 31488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 58983, 0, 3,
                                                                       58443, 31158, 58488,
                                                                       14910, 14946, 31638,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59073, 0, 3,
                                                                       58488, 31188, 58533,
                                                                       14946, 14982, 31698,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59163, 0, 3,
                                                                       58533, 31218, 58578,
                                                                       14982, 15018, 31758,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59253, 0, 3,
                                                                       58578, 31248, 58623,
                                                                       15018, 15054, 31818,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59343, 0, 3,
                                                                       58623, 31278, 58668,
                                                                       15054, 15090, 31878,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59433, 0, 3,
                                                                       58668, 31308, 58713,
                                                                       15090, 15126, 31938,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59523, 0, 3,
                                                                       58713, 31338, 58758,
                                                                       15126, 15162, 31998,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59613, 0, 3,
                                                                       58758, 31368, 58803,
                                                                       15162, 15198, 32058,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59703, 0, 3,
                                                                       58803, 31398, 58848,
                                                                       15198, 15234, 32118,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59793, 0, 3,
                                                                       58848, 31428, 58893,
                                                                       15234, 15270, 32178,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59883, 0, 3,
                                                                       58893, 31458, 58938,
                                                                       15270, 15306, 32238,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 59973, 0, 3,
                                                                       58983, 31638, 59073,
                                                                       15378, 15438, 32498,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60123, 0, 3,
                                                                       59073, 31698, 59163,
                                                                       15438, 15498, 32598,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60273, 0, 3,
                                                                       59163, 31758, 59253,
                                                                       15498, 15558, 32698,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60423, 0, 3,
                                                                       59253, 31818, 59343,
                                                                       15558, 15618, 32798,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60573, 0, 3,
                                                                       59343, 31878, 59433,
                                                                       15618, 15678, 32898,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60723, 0, 3,
                                                                       59433, 31938, 59523,
                                                                       15678, 15738, 32998,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60873, 0, 3,
                                                                       59523, 31998, 59613,
                                                                       15738, 15798, 33098,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61023, 0, 3,
                                                                       59613, 32058, 59703,
                                                                       15798, 15858, 33198,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61173, 0, 3,
                                                                       59703, 32118, 59793,
                                                                       15858, 15918, 33298,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61323, 0, 3,
                                                                       59793, 32178, 59883,
                                                                       15918, 15978, 33398,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 61473, 0, 3,
                                                                       59973, 32498, 60123,
                                                                       16098, 16188, 33798,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 61698, 0, 3,
                                                                       60123, 32598, 60273,
                                                                       16188, 16278, 33948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 61923, 0, 3,
                                                                       60273, 32698, 60423,
                                                                       16278, 16368, 34098,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62148, 0, 3,
                                                                       60423, 32798, 60573,
                                                                       16368, 16458, 34248,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62373, 0, 3,
                                                                       60573, 32898, 60723,
                                                                       16458, 16548, 34398,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62598, 0, 3,
                                                                       60723, 32998, 60873,
                                                                       16548, 16638, 34548,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62823, 0, 3,
                                                                       60873, 33098, 61023,
                                                                       16638, 16728, 34698,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63048, 0, 3,
                                                                       61023, 33198, 61173,
                                                                       16728, 16818, 34848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63273, 0, 3,
                                                                       61173, 33298, 61323,
                                                                       16818, 16908, 34998,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 63498, 0, 3,
                                                                       61473, 33798, 61698,
                                                                       17088, 17214, 35568,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 63813, 0, 3,
                                                                       61698, 33948, 61923,
                                                                       17214, 17340, 35778,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 64128, 0, 3,
                                                                       61923, 34098, 62148,
                                                                       17340, 17466, 35988,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 64443, 0, 3,
                                                                       62148, 34248, 62373,
                                                                       17466, 17592, 36198,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 64758, 0, 3,
                                                                       62373, 34398, 62598,
                                                                       17592, 17718, 36408,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65073, 0, 3,
                                                                       62598, 34548, 62823,
                                                                       17718, 17844, 36618,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65388, 0, 3,
                                                                       62823, 34698, 63048,
                                                                       17844, 17970, 36828,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65703, 0, 3,
                                                                       63048, 34848, 63273,
                                                                       17970, 18096, 37038,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 66018, 0, 3,
                                                                       63498, 35568, 63813,
                                                                       18348, 18516, 37808,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 66438, 0, 3,
                                                                       63813, 35778, 64128,
                                                                       18516, 18684, 38088,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 66858, 0, 3,
                                                                       64128, 35988, 64443,
                                                                       18684, 18852, 38368,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 67278, 0, 3,
                                                                       64443, 36198, 64758,
                                                                       18852, 19020, 38648,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 67698, 0, 3,
                                                                       64758, 36408, 65073,
                                                                       19020, 19188, 38928,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 68118, 0, 3,
                                                                       65073, 36618, 65388,
                                                                       19188, 19356, 39208,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 68538, 0, 3,
                                                                       65388, 36828, 65703,
                                                                       19356, 19524, 39488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 68958, 0, 3,
                                                                       66018, 37808, 66438,
                                                                       19860, 20076, 40488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 69498, 0, 3,
                                                                       66438, 38088, 66858,
                                                                       20076, 20292, 40848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 70038, 0, 3,
                                                                       66858, 38368, 67278,
                                                                       20292, 20508, 41208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 70578, 0, 3,
                                                                       67278, 38648, 67698,
                                                                       20508, 20724, 41568,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 71118, 0, 3,
                                                                       67698, 38928, 68118,
                                                                       20724, 20940, 41928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 71658, 0, 3,
                                                                       68118, 39208, 68538,
                                                                       20940, 21156, 42288,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 72198, 0, 3,
                                                                       68958, 40488, 69498,
                                                                       21588, 21858, 43548,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 72873, 0, 3,
                                                                       69498, 40848, 70038,
                                                                       21858, 22128, 43998,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 73548, 0, 3,
                                                                       70038, 41208, 70578,
                                                                       22128, 22398, 44448,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 74223, 0, 3,
                                                                       70578, 41568, 71118,
                                                                       22398, 22668, 44898,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 74898, 0, 3,
                                                                       71118, 41928, 71658,
                                                                       22668, 22938, 45348,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 75573, 0, 3,
                                                                       72198, 43548, 72873,
                                                                       23478, 23808, 46898,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 76398, 0, 3,
                                                                       72873, 43998, 73548,
                                                                       23808, 24138, 47448,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 77223, 0, 3,
                                                                       73548, 44448, 74223,
                                                                       24138, 24468, 47998,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 78048, 0, 3,
                                                                       74223, 44898, 74898,
                                                                       24468, 24798, 48548,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 78873, 0, 3,
                                                                       75573, 46898, 76398,
                                                                       25458, 25854, 50418,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 79863, 0, 3,
                                                                       76398, 47448, 77223,
                                                                       25854, 26250, 51078,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 80853, 0, 3,
                                                                       77223, 47998, 78048,
                                                                       26250, 26646, 51738,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 81843, 0, 3,
                                                                       78873, 50418, 79863,
                                                                       27438, 27906, 53958,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 83013, 0, 3,
                                                                       79863, 51078, 80853,
                                                                       27906, 28374, 54738,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 84183, 0, 3,
                                                                       81843, 53958, 83013,
                                                                       29310, 29856, 57338,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85548, 3, 30948,
                                                                       30958, 58248, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85569, 3, 30958,
                                                                       30968, 58263, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85590, 3, 30968,
                                                                       30978, 58278, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85611, 3, 30978,
                                                                       30988, 58293, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85632, 3, 30988,
                                                                       30998, 58308, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85653, 3, 30998,
                                                                       31008, 58323, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85674, 3, 31008,
                                                                       31018, 58338, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85695, 3, 31018,
                                                                       31028, 58353, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85716, 3, 31028,
                                                                       31038, 58368, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85737, 3, 31038,
                                                                       31048, 58383, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85758, 3, 31048,
                                                                       31058, 58398, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85779, 3, 31058,
                                                                       31068, 58413, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85800, 3, 31068,
                                                                       31078, 58428, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 85821, 0, 3,
                                                                       85548, 58248, 85569,
                                                                       31098, 31128, 58443,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 85884, 0, 3,
                                                                       85569, 58263, 85590,
                                                                       31128, 31158, 58488,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 85947, 0, 3,
                                                                       85590, 58278, 85611,
                                                                       31158, 31188, 58533,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86010, 0, 3,
                                                                       85611, 58293, 85632,
                                                                       31188, 31218, 58578,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86073, 0, 3,
                                                                       85632, 58308, 85653,
                                                                       31218, 31248, 58623,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86136, 0, 3,
                                                                       85653, 58323, 85674,
                                                                       31248, 31278, 58668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86199, 0, 3,
                                                                       85674, 58338, 85695,
                                                                       31278, 31308, 58713,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86262, 0, 3,
                                                                       85695, 58353, 85716,
                                                                       31308, 31338, 58758,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86325, 0, 3,
                                                                       85716, 58368, 85737,
                                                                       31338, 31368, 58803,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86388, 0, 3,
                                                                       85737, 58383, 85758,
                                                                       31368, 31398, 58848,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86451, 0, 3,
                                                                       85758, 58398, 85779,
                                                                       31398, 31428, 58893,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86514, 0, 3,
                                                                       85779, 58413, 85800,
                                                                       31428, 31458, 58938,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 86577, 0, 3,
                                                                       85821, 58443, 85884,
                                                                       31518, 31578, 58983,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 86703, 0, 3,
                                                                       85884, 58488, 85947,
                                                                       31578, 31638, 59073,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 86829, 0, 3,
                                                                       85947, 58533, 86010,
                                                                       31638, 31698, 59163,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 86955, 0, 3,
                                                                       86010, 58578, 86073,
                                                                       31698, 31758, 59253,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87081, 0, 3,
                                                                       86073, 58623, 86136,
                                                                       31758, 31818, 59343,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87207, 0, 3,
                                                                       86136, 58668, 86199,
                                                                       31818, 31878, 59433,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87333, 0, 3,
                                                                       86199, 58713, 86262,
                                                                       31878, 31938, 59523,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87459, 0, 3,
                                                                       86262, 58758, 86325,
                                                                       31938, 31998, 59613,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87585, 0, 3,
                                                                       86325, 58803, 86388,
                                                                       31998, 32058, 59703,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87711, 0, 3,
                                                                       86388, 58848, 86451,
                                                                       32058, 32118, 59793,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87837, 0, 3,
                                                                       86451, 58893, 86514,
                                                                       32118, 32178, 59883,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 87963, 0, 3,
                                                                       86577, 58983, 86703,
                                                                       32298, 32398, 59973,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 88173, 0, 3,
                                                                       86703, 59073, 86829,
                                                                       32398, 32498, 60123,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 88383, 0, 3,
                                                                       86829, 59163, 86955,
                                                                       32498, 32598, 60273,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 88593, 0, 3,
                                                                       86955, 59253, 87081,
                                                                       32598, 32698, 60423,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 88803, 0, 3,
                                                                       87081, 59343, 87207,
                                                                       32698, 32798, 60573,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 89013, 0, 3,
                                                                       87207, 59433, 87333,
                                                                       32798, 32898, 60723,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 89223, 0, 3,
                                                                       87333, 59523, 87459,
                                                                       32898, 32998, 60873,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 89433, 0, 3,
                                                                       87459, 59613, 87585,
                                                                       32998, 33098, 61023,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 89643, 0, 3,
                                                                       87585, 59703, 87711,
                                                                       33098, 33198, 61173,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 89853, 0, 3,
                                                                       87711, 59793, 87837,
                                                                       33198, 33298, 61323,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 90063, 0, 3,
                                                                       87963, 59973, 88173,
                                                                       33498, 33648, 61473,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 90378, 0, 3,
                                                                       88173, 60123, 88383,
                                                                       33648, 33798, 61698,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 90693, 0, 3,
                                                                       88383, 60273, 88593,
                                                                       33798, 33948, 61923,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 91008, 0, 3,
                                                                       88593, 60423, 88803,
                                                                       33948, 34098, 62148,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 91323, 0, 3,
                                                                       88803, 60573, 89013,
                                                                       34098, 34248, 62373,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 91638, 0, 3,
                                                                       89013, 60723, 89223,
                                                                       34248, 34398, 62598,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 91953, 0, 3,
                                                                       89223, 60873, 89433,
                                                                       34398, 34548, 62823,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 92268, 0, 3,
                                                                       89433, 61023, 89643,
                                                                       34548, 34698, 63048,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 92583, 0, 3,
                                                                       89643, 61173, 89853,
                                                                       34698, 34848, 63273,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 92898, 0, 3,
                                                                       90063, 61473, 90378,
                                                                       35148, 35358, 63498,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 93339, 0, 3,
                                                                       90378, 61698, 90693,
                                                                       35358, 35568, 63813,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 93780, 0, 3,
                                                                       90693, 61923, 91008,
                                                                       35568, 35778, 64128,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 94221, 0, 3,
                                                                       91008, 62148, 91323,
                                                                       35778, 35988, 64443,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 94662, 0, 3,
                                                                       91323, 62373, 91638,
                                                                       35988, 36198, 64758,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 95103, 0, 3,
                                                                       91638, 62598, 91953,
                                                                       36198, 36408, 65073,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 95544, 0, 3,
                                                                       91953, 62823, 92268,
                                                                       36408, 36618, 65388,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 95985, 0, 3,
                                                                       92268, 63048, 92583,
                                                                       36618, 36828, 65703,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 96426, 0, 3,
                                                                       92898, 63498, 93339,
                                                                       37248, 37528, 66018,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 97014, 0, 3,
                                                                       93339, 63813, 93780,
                                                                       37528, 37808, 66438,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 97602, 0, 3,
                                                                       93780, 64128, 94221,
                                                                       37808, 38088, 66858,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 98190, 0, 3,
                                                                       94221, 64443, 94662,
                                                                       38088, 38368, 67278,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 98778, 0, 3,
                                                                       94662, 64758, 95103,
                                                                       38368, 38648, 67698,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 99366, 0, 3,
                                                                       95103, 65073, 95544,
                                                                       38648, 38928, 68118,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 99954, 0, 3,
                                                                       95544, 65388, 95985,
                                                                       38928, 39208, 68538,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 100542, 0, 3,
                                                                       96426, 66018, 97014,
                                                                       39768, 40128, 68958,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 101298, 0, 3,
                                                                       97014, 66438, 97602,
                                                                       40128, 40488, 69498,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 102054, 0, 3,
                                                                       97602, 66858, 98190,
                                                                       40488, 40848, 70038,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 102810, 0, 3,
                                                                       98190, 67278, 98778,
                                                                       40848, 41208, 70578,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 103566, 0, 3,
                                                                       98778, 67698, 99366,
                                                                       41208, 41568, 71118,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 104322, 0, 3,
                                                                       99366, 68118, 99954,
                                                                       41568, 41928, 71658,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 105078, 0, 3,
                                                                       100542, 68958, 101298,
                                                                       42648, 43098, 72198,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 106023, 0, 3,
                                                                       101298, 69498, 102054,
                                                                       43098, 43548, 72873,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 106968, 0, 3,
                                                                       102054, 70038, 102810,
                                                                       43548, 43998, 73548,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 107913, 0, 3,
                                                                       102810, 70578, 103566,
                                                                       43998, 44448, 74223,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 108858, 0, 3,
                                                                       103566, 71118, 104322,
                                                                       44448, 44898, 74898,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 109803, 0, 3,
                                                                       105078, 72198, 106023,
                                                                       45798, 46348, 75573,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 110958, 0, 3,
                                                                       106023, 72873, 106968,
                                                                       46348, 46898, 76398,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 112113, 0, 3,
                                                                       106968, 73548, 107913,
                                                                       46898, 47448, 77223,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 113268, 0, 3,
                                                                       107913, 74223, 108858,
                                                                       47448, 47998, 78048,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 114423, 0, 3,
                                                                       109803, 75573, 110958,
                                                                       49098, 49758, 78873,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 115809, 0, 3,
                                                                       110958, 76398, 112113,
                                                                       49758, 50418, 79863,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 117195, 0, 3,
                                                                       112113, 77223, 113268,
                                                                       50418, 51078, 80853,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 118581, 0, 3,
                                                                       114423, 78873, 115809,
                                                                       52398, 53178, 81843,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 120219, 0, 3,
                                                                       115809, 79863, 117195,
                                                                       53178, 53958, 83013,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 121857, 0, 3,
                                                                       118581, 81843, 120219,
                                                                       55518, 56428, 84183,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 123768, 96426, 588, ncols);

                    simdfunc::contract_primitives(buffer, 124664, 100542, 756, ncols);

                    simdfunc::contract_primitives(buffer, 125816, 105078, 945, ncols);

                    simdfunc::contract_primitives(buffer, 127256, 109803, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 129016, 114423, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 131128, 118581, 1638, ncols);

                    simdfunc::contract_primitives(buffer, 133624, 121857, 1911, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 124356, 123768, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 125420, 124664, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 126761, 125816, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 128411, 127256, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 130402, 129016, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 132766, 131128, 78, 1, nmax);

        simdtrf::transform_h_inner(buffer, 135535, 133624, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 136536, 124356, 125420, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 137460, 125420, 126761, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 138648, 126761, 128411, 11,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 140133, 128411, 130402, 11,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 141948, 130402, 132766, 11,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 144126, 132766, 135535, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 146700, 136536, 137460, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 148548, 137460, 138648, 11,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 150924, 138648, 140133, 11,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 153894, 140133, 141948, 11,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 157524, 141948, 144126, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 161880, 146700, 148548, 11,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 164960, 148548, 150924, 11,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 168920, 150924, 153894, 11,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 173870, 153894, 157524, 11,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 179920, 161880, 164960, 11,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 184540, 164960, 168920, 11,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 190480, 168920, 173870, 11,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 197905, 179920, 184540, 11,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 204373, 184540, 190480, 11,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 212689, 197905, 204373, 11,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 221313, 212689, 28, 11, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 221313, 143, nmax);
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
