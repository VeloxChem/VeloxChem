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


#include "SimdThreeCenterElectronRepulsionRecIIG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iig_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iig_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 155933, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1521 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 155933, 73718, 8757, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 16,
                                                             ncols, fj, 6, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3648, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3651, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3654, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3657, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3660, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3663, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3666, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3669, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3672, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3675, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3678, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3681, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3684, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3687, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3690, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3693, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3702, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3711, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3720, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3729, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3738, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3747, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3756, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3765, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3774, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3783, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3792, 3, 21, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3801, 3, 22, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3810, 3, 23, 70,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3819, 3, 31, 85,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3837, 3, 34, 91,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3855, 3, 37, 97,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3873, 3, 40, 103,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3891, 3, 43, 109,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3909, 3, 46, 115,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3927, 3, 49, 121,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3945, 3, 52, 127,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3963, 3, 55, 133,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3981, 3, 58, 139,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3999, 3, 61, 145,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4017, 3, 64, 151,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4035, 3, 67, 157,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4053, 3, 85, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4083, 3, 91, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4113, 3, 97, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4143, 3, 103, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4173, 3, 109, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4203, 3, 115, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4233, 3, 121, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4263, 3, 127, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4293, 3, 133, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4323, 3, 139, 273,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4353, 3, 145, 283,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4383, 3, 151, 293,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4413, 3, 183, 333,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4458, 3, 193, 348,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4503, 3, 203, 363,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4548, 3, 213, 378,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4593, 3, 223, 393,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4638, 3, 233, 408,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4683, 3, 243, 423,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4728, 3, 253, 438,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4773, 3, 263, 453,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4818, 3, 273, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4863, 3, 283, 483,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4908, 3, 333, 540,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4971, 3, 348, 561,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5034, 3, 363, 582,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5097, 3, 378, 603,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5160, 3, 393, 624,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5223, 3, 408, 645,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5286, 3, 423, 666,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5349, 3, 438, 687,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5412, 3, 453, 708,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5475, 3, 468, 729,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5538, 3, 540, 806,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5622, 3, 561, 834,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5706, 3, 582, 862,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5790, 3, 603, 890,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5874, 3, 624, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5958, 3, 645, 946,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6042, 3, 666, 974,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6126, 3, 687,
                                                                       1002, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6210, 3, 708,
                                                                       1030, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6294, 3, 806,
                                                                       1130, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6402, 3, 834,
                                                                       1166, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6510, 3, 862,
                                                                       1202, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6618, 3, 890,
                                                                       1238, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6726, 3, 918,
                                                                       1274, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6834, 3, 946,
                                                                       1310, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6942, 3, 974,
                                                                       1346, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7050, 3, 1002,
                                                                       1382, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7158, 3, 1130,
                                                                       1508, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7293, 3, 1166,
                                                                       1553, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7428, 3, 1202,
                                                                       1598, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7563, 3, 1238,
                                                                       1643, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7698, 3, 1274,
                                                                       1688, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7833, 3, 1310,
                                                                       1733, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7968, 3, 1346,
                                                                       1778, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8103, 3, 1508,
                                                                       1933, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8268, 3, 1553,
                                                                       1988, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8433, 3, 1598,
                                                                       2043, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8598, 3, 1643,
                                                                       2098, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8763, 3, 1688,
                                                                       2153, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8928, 3, 1733,
                                                                       2208, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9093, 3, 1933,
                                                                       2395, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9291, 3, 1988,
                                                                       2461, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9489, 3, 2043,
                                                                       2527, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9687, 3, 2098,
                                                                       2593, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9885, 3, 2153,
                                                                       2659, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10083, 3, 2395,
                                                                       2881, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10317, 3, 2461,
                                                                       2959, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10551, 3, 2527,
                                                                       3037, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10785, 3, 2593,
                                                                       3115, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 11019, 3, 2881,
                                                                       3375, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 11292, 3, 2959,
                                                                       3466, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 11565, 3, 3037,
                                                                       3557, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11838, 3, 8, 9,
                                                                       3648, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11844, 3, 9, 10,
                                                                       3651, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11850, 3, 10, 11,
                                                                       3654, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11856, 3, 11, 12,
                                                                       3657, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11862, 3, 12, 13,
                                                                       3660, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11868, 3, 13, 14,
                                                                       3663, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11874, 3, 14, 15,
                                                                       3666, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11880, 3, 15, 16,
                                                                       3669, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11886, 3, 16, 17,
                                                                       3672, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11892, 3, 17, 18,
                                                                       3675, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11898, 3, 18, 19,
                                                                       3678, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11904, 3, 19, 20,
                                                                       3681, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11910, 3, 20, 21,
                                                                       3684, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11916, 3, 21, 22,
                                                                       3687, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11922, 3, 22, 23,
                                                                       3690, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11928, 0, 3,
                                                                       11838, 3648, 11844, 3693,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11946, 0, 3,
                                                                       11844, 3651, 11850, 3702,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11964, 0, 3,
                                                                       11850, 3654, 11856, 3711,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11982, 0, 3,
                                                                       11856, 3657, 11862, 3720,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12000, 0, 3,
                                                                       11862, 3660, 11868, 3729,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12018, 0, 3,
                                                                       11868, 3663, 11874, 3738,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12036, 0, 3,
                                                                       11874, 3666, 11880, 3747,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12054, 0, 3,
                                                                       11880, 3669, 11886, 3756,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12072, 0, 3,
                                                                       11886, 3672, 11892, 3765,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12090, 0, 3,
                                                                       11892, 3675, 11898, 3774,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12108, 0, 3,
                                                                       11898, 3678, 11904, 3783,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12126, 0, 3,
                                                                       11904, 3681, 11910, 3792,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12144, 0, 3,
                                                                       11910, 3684, 11916, 3801,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12162, 0, 3,
                                                                       11916, 3687, 11922, 3810,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12180, 0, 3,
                                                                       11928, 3693, 11946, 73,
                                                                       79, 3819, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12216, 0, 3,
                                                                       11946, 3702, 11964, 79,
                                                                       85, 3837, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12252, 0, 3,
                                                                       11964, 3711, 11982, 85,
                                                                       91, 3855, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12288, 0, 3,
                                                                       11982, 3720, 12000, 91,
                                                                       97, 3873, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12324, 0, 3,
                                                                       12000, 3729, 12018, 97,
                                                                       103, 3891, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12360, 0, 3,
                                                                       12018, 3738, 12036, 103,
                                                                       109, 3909, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12396, 0, 3,
                                                                       12036, 3747, 12054, 109,
                                                                       115, 3927, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12432, 0, 3,
                                                                       12054, 3756, 12072, 115,
                                                                       121, 3945, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12468, 0, 3,
                                                                       12072, 3765, 12090, 121,
                                                                       127, 3963, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12504, 0, 3,
                                                                       12090, 3774, 12108, 127,
                                                                       133, 3981, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12540, 0, 3,
                                                                       12108, 3783, 12126, 133,
                                                                       139, 3999, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12576, 0, 3,
                                                                       12126, 3792, 12144, 139,
                                                                       145, 4017, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12612, 0, 3,
                                                                       12144, 3801, 12162, 145,
                                                                       151, 4035, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12648, 0, 3,
                                                                       12180, 3819, 12216, 163,
                                                                       173, 4053, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12708, 0, 3,
                                                                       12216, 3837, 12252, 173,
                                                                       183, 4083, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12768, 0, 3,
                                                                       12252, 3855, 12288, 183,
                                                                       193, 4113, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12828, 0, 3,
                                                                       12288, 3873, 12324, 193,
                                                                       203, 4143, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12888, 0, 3,
                                                                       12324, 3891, 12360, 203,
                                                                       213, 4173, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12948, 0, 3,
                                                                       12360, 3909, 12396, 213,
                                                                       223, 4203, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13008, 0, 3,
                                                                       12396, 3927, 12432, 223,
                                                                       233, 4233, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13068, 0, 3,
                                                                       12432, 3945, 12468, 233,
                                                                       243, 4263, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13128, 0, 3,
                                                                       12468, 3963, 12504, 243,
                                                                       253, 4293, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13188, 0, 3,
                                                                       12504, 3981, 12540, 253,
                                                                       263, 4323, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13248, 0, 3,
                                                                       12540, 3999, 12576, 263,
                                                                       273, 4353, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13308, 0, 3,
                                                                       12576, 4017, 12612, 273,
                                                                       283, 4383, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13368, 0, 3,
                                                                       12648, 4053, 12708, 303,
                                                                       318, 4413, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13458, 0, 3,
                                                                       12708, 4083, 12768, 318,
                                                                       333, 4458, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13548, 0, 3,
                                                                       12768, 4113, 12828, 333,
                                                                       348, 4503, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13638, 0, 3,
                                                                       12828, 4143, 12888, 348,
                                                                       363, 4548, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13728, 0, 3,
                                                                       12888, 4173, 12948, 363,
                                                                       378, 4593, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13818, 0, 3,
                                                                       12948, 4203, 13008, 378,
                                                                       393, 4638, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13908, 0, 3,
                                                                       13008, 4233, 13068, 393,
                                                                       408, 4683, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13998, 0, 3,
                                                                       13068, 4263, 13128, 408,
                                                                       423, 4728, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14088, 0, 3,
                                                                       13128, 4293, 13188, 423,
                                                                       438, 4773, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14178, 0, 3,
                                                                       13188, 4323, 13248, 438,
                                                                       453, 4818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14268, 0, 3,
                                                                       13248, 4353, 13308, 453,
                                                                       468, 4863, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14358, 0, 3,
                                                                       13368, 4413, 13458, 498,
                                                                       519, 4908, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14484, 0, 3,
                                                                       13458, 4458, 13548, 519,
                                                                       540, 4971, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14610, 0, 3,
                                                                       13548, 4503, 13638, 540,
                                                                       561, 5034, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14736, 0, 3,
                                                                       13638, 4548, 13728, 561,
                                                                       582, 5097, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14862, 0, 3,
                                                                       13728, 4593, 13818, 582,
                                                                       603, 5160, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14988, 0, 3,
                                                                       13818, 4638, 13908, 603,
                                                                       624, 5223, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15114, 0, 3,
                                                                       13908, 4683, 13998, 624,
                                                                       645, 5286, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15240, 0, 3,
                                                                       13998, 4728, 14088, 645,
                                                                       666, 5349, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15366, 0, 3,
                                                                       14088, 4773, 14178, 666,
                                                                       687, 5412, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15492, 0, 3,
                                                                       14178, 4818, 14268, 687,
                                                                       708, 5475, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15618, 0, 3,
                                                                       14358, 4908, 14484, 750,
                                                                       778, 5538, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15786, 0, 3,
                                                                       14484, 4971, 14610, 778,
                                                                       806, 5622, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15954, 0, 3,
                                                                       14610, 5034, 14736, 806,
                                                                       834, 5706, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16122, 0, 3,
                                                                       14736, 5097, 14862, 834,
                                                                       862, 5790, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16290, 0, 3,
                                                                       14862, 5160, 14988, 862,
                                                                       890, 5874, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16458, 0, 3,
                                                                       14988, 5223, 15114, 890,
                                                                       918, 5958, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16626, 0, 3,
                                                                       15114, 5286, 15240, 918,
                                                                       946, 6042, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16794, 0, 3,
                                                                       15240, 5349, 15366, 946,
                                                                       974, 6126, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16962, 0, 3,
                                                                       15366, 5412, 15492, 974,
                                                                       1002, 6210, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17130, 0, 3,
                                                                       15618, 5538, 15786, 1058,
                                                                       1094, 6294, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17346, 0, 3,
                                                                       15786, 5622, 15954, 1094,
                                                                       1130, 6402, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17562, 0, 3,
                                                                       15954, 5706, 16122, 1130,
                                                                       1166, 6510, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17778, 0, 3,
                                                                       16122, 5790, 16290, 1166,
                                                                       1202, 6618, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17994, 0, 3,
                                                                       16290, 5874, 16458, 1202,
                                                                       1238, 6726, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18210, 0, 3,
                                                                       16458, 5958, 16626, 1238,
                                                                       1274, 6834, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18426, 0, 3,
                                                                       16626, 6042, 16794, 1274,
                                                                       1310, 6942, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18642, 0, 3,
                                                                       16794, 6126, 16962, 1310,
                                                                       1346, 7050, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18858, 0, 3,
                                                                       17130, 6294, 17346, 1418,
                                                                       1463, 7158, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19128, 0, 3,
                                                                       17346, 6402, 17562, 1463,
                                                                       1508, 7293, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19398, 0, 3,
                                                                       17562, 6510, 17778, 1508,
                                                                       1553, 7428, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19668, 0, 3,
                                                                       17778, 6618, 17994, 1553,
                                                                       1598, 7563, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19938, 0, 3,
                                                                       17994, 6726, 18210, 1598,
                                                                       1643, 7698, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20208, 0, 3,
                                                                       18210, 6834, 18426, 1643,
                                                                       1688, 7833, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20478, 0, 3,
                                                                       18426, 6942, 18642, 1688,
                                                                       1733, 7968, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20748, 0, 3,
                                                                       18858, 7158, 19128, 1823,
                                                                       1878, 8103, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 21078, 0, 3,
                                                                       19128, 7293, 19398, 1878,
                                                                       1933, 8268, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 21408, 0, 3,
                                                                       19398, 7428, 19668, 1933,
                                                                       1988, 8433, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 21738, 0, 3,
                                                                       19668, 7563, 19938, 1988,
                                                                       2043, 8598, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22068, 0, 3,
                                                                       19938, 7698, 20208, 2043,
                                                                       2098, 8763, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22398, 0, 3,
                                                                       20208, 7833, 20478, 2098,
                                                                       2153, 8928, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 22728, 0, 3,
                                                                       20748, 8103, 21078, 2263,
                                                                       2329, 9093, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 23124, 0, 3,
                                                                       21078, 8268, 21408, 2329,
                                                                       2395, 9291, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 23520, 0, 3,
                                                                       21408, 8433, 21738, 2395,
                                                                       2461, 9489, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 23916, 0, 3,
                                                                       21738, 8598, 22068, 2461,
                                                                       2527, 9687, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 24312, 0, 3,
                                                                       22068, 8763, 22398, 2527,
                                                                       2593, 9885, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 24708, 0, 3,
                                                                       22728, 9093, 23124, 2725,
                                                                       2803, 10083, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 25176, 0, 3,
                                                                       23124, 9291, 23520, 2803,
                                                                       2881, 10317, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 25644, 0, 3,
                                                                       23520, 9489, 23916, 2881,
                                                                       2959, 10551, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 26112, 0, 3,
                                                                       23916, 9687, 24312, 2959,
                                                                       3037, 10785, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 26580, 0, 3,
                                                                       24708, 10083, 25176, 3193,
                                                                       3284, 11019, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 27126, 0, 3,
                                                                       25176, 10317, 25644, 3284,
                                                                       3375, 11292, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 27672, 0, 3,
                                                                       25644, 10551, 26112, 3375,
                                                                       3466, 11565, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28218, 3, 3648,
                                                                       3651, 11850, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28228, 3, 3651,
                                                                       3654, 11856, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28238, 3, 3654,
                                                                       3657, 11862, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28248, 3, 3657,
                                                                       3660, 11868, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28258, 3, 3660,
                                                                       3663, 11874, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28268, 3, 3663,
                                                                       3666, 11880, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28278, 3, 3666,
                                                                       3669, 11886, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28288, 3, 3669,
                                                                       3672, 11892, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28298, 3, 3672,
                                                                       3675, 11898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28308, 3, 3675,
                                                                       3678, 11904, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28318, 3, 3678,
                                                                       3681, 11910, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28328, 3, 3681,
                                                                       3684, 11916, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28338, 3, 3684,
                                                                       3687, 11922, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28348, 0, 3,
                                                                       28218, 11850, 28228, 3693,
                                                                       3702, 11964, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28378, 0, 3,
                                                                       28228, 11856, 28238, 3702,
                                                                       3711, 11982, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28408, 0, 3,
                                                                       28238, 11862, 28248, 3711,
                                                                       3720, 12000, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28438, 0, 3,
                                                                       28248, 11868, 28258, 3720,
                                                                       3729, 12018, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28468, 0, 3,
                                                                       28258, 11874, 28268, 3729,
                                                                       3738, 12036, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28498, 0, 3,
                                                                       28268, 11880, 28278, 3738,
                                                                       3747, 12054, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28528, 0, 3,
                                                                       28278, 11886, 28288, 3747,
                                                                       3756, 12072, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28558, 0, 3,
                                                                       28288, 11892, 28298, 3756,
                                                                       3765, 12090, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28588, 0, 3,
                                                                       28298, 11898, 28308, 3765,
                                                                       3774, 12108, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28618, 0, 3,
                                                                       28308, 11904, 28318, 3774,
                                                                       3783, 12126, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28648, 0, 3,
                                                                       28318, 11910, 28328, 3783,
                                                                       3792, 12144, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28678, 0, 3,
                                                                       28328, 11916, 28338, 3792,
                                                                       3801, 12162, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28708, 0, 3,
                                                                       28348, 11964, 28378, 3819,
                                                                       3837, 12252, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28768, 0, 3,
                                                                       28378, 11982, 28408, 3837,
                                                                       3855, 12288, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28828, 0, 3,
                                                                       28408, 12000, 28438, 3855,
                                                                       3873, 12324, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28888, 0, 3,
                                                                       28438, 12018, 28468, 3873,
                                                                       3891, 12360, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28948, 0, 3,
                                                                       28468, 12036, 28498, 3891,
                                                                       3909, 12396, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29008, 0, 3,
                                                                       28498, 12054, 28528, 3909,
                                                                       3927, 12432, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29068, 0, 3,
                                                                       28528, 12072, 28558, 3927,
                                                                       3945, 12468, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29128, 0, 3,
                                                                       28558, 12090, 28588, 3945,
                                                                       3963, 12504, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29188, 0, 3,
                                                                       28588, 12108, 28618, 3963,
                                                                       3981, 12540, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29248, 0, 3,
                                                                       28618, 12126, 28648, 3981,
                                                                       3999, 12576, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29308, 0, 3,
                                                                       28648, 12144, 28678, 3999,
                                                                       4017, 12612, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29368, 0, 3,
                                                                       28708, 12252, 28768, 4053,
                                                                       4083, 12768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29468, 0, 3,
                                                                       28768, 12288, 28828, 4083,
                                                                       4113, 12828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29568, 0, 3,
                                                                       28828, 12324, 28888, 4113,
                                                                       4143, 12888, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29668, 0, 3,
                                                                       28888, 12360, 28948, 4143,
                                                                       4173, 12948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29768, 0, 3,
                                                                       28948, 12396, 29008, 4173,
                                                                       4203, 13008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29868, 0, 3,
                                                                       29008, 12432, 29068, 4203,
                                                                       4233, 13068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29968, 0, 3,
                                                                       29068, 12468, 29128, 4233,
                                                                       4263, 13128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30068, 0, 3,
                                                                       29128, 12504, 29188, 4263,
                                                                       4293, 13188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30168, 0, 3,
                                                                       29188, 12540, 29248, 4293,
                                                                       4323, 13248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30268, 0, 3,
                                                                       29248, 12576, 29308, 4323,
                                                                       4353, 13308, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30368, 0, 3,
                                                                       29368, 12768, 29468, 4413,
                                                                       4458, 13548, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30518, 0, 3,
                                                                       29468, 12828, 29568, 4458,
                                                                       4503, 13638, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30668, 0, 3,
                                                                       29568, 12888, 29668, 4503,
                                                                       4548, 13728, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30818, 0, 3,
                                                                       29668, 12948, 29768, 4548,
                                                                       4593, 13818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30968, 0, 3,
                                                                       29768, 13008, 29868, 4593,
                                                                       4638, 13908, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31118, 0, 3,
                                                                       29868, 13068, 29968, 4638,
                                                                       4683, 13998, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31268, 0, 3,
                                                                       29968, 13128, 30068, 4683,
                                                                       4728, 14088, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31418, 0, 3,
                                                                       30068, 13188, 30168, 4728,
                                                                       4773, 14178, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31568, 0, 3,
                                                                       30168, 13248, 30268, 4773,
                                                                       4818, 14268, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 31718, 0, 3,
                                                                       30368, 13548, 30518, 4908,
                                                                       4971, 14610, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 31928, 0, 3,
                                                                       30518, 13638, 30668, 4971,
                                                                       5034, 14736, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32138, 0, 3,
                                                                       30668, 13728, 30818, 5034,
                                                                       5097, 14862, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32348, 0, 3,
                                                                       30818, 13818, 30968, 5097,
                                                                       5160, 14988, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32558, 0, 3,
                                                                       30968, 13908, 31118, 5160,
                                                                       5223, 15114, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32768, 0, 3,
                                                                       31118, 13998, 31268, 5223,
                                                                       5286, 15240, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32978, 0, 3,
                                                                       31268, 14088, 31418, 5286,
                                                                       5349, 15366, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33188, 0, 3,
                                                                       31418, 14178, 31568, 5349,
                                                                       5412, 15492, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 33398, 0, 3,
                                                                       31718, 14610, 31928, 5538,
                                                                       5622, 15954, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 33678, 0, 3,
                                                                       31928, 14736, 32138, 5622,
                                                                       5706, 16122, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 33958, 0, 3,
                                                                       32138, 14862, 32348, 5706,
                                                                       5790, 16290, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34238, 0, 3,
                                                                       32348, 14988, 32558, 5790,
                                                                       5874, 16458, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34518, 0, 3,
                                                                       32558, 15114, 32768, 5874,
                                                                       5958, 16626, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34798, 0, 3,
                                                                       32768, 15240, 32978, 5958,
                                                                       6042, 16794, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35078, 0, 3,
                                                                       32978, 15366, 33188, 6042,
                                                                       6126, 16962, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 35358, 0, 3,
                                                                       33398, 15954, 33678, 6294,
                                                                       6402, 17562, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 35718, 0, 3,
                                                                       33678, 16122, 33958, 6402,
                                                                       6510, 17778, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36078, 0, 3,
                                                                       33958, 16290, 34238, 6510,
                                                                       6618, 17994, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36438, 0, 3,
                                                                       34238, 16458, 34518, 6618,
                                                                       6726, 18210, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36798, 0, 3,
                                                                       34518, 16626, 34798, 6726,
                                                                       6834, 18426, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 37158, 0, 3,
                                                                       34798, 16794, 35078, 6834,
                                                                       6942, 18642, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 37518, 0, 3,
                                                                       35358, 17562, 35718, 7158,
                                                                       7293, 19398, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 37968, 0, 3,
                                                                       35718, 17778, 36078, 7293,
                                                                       7428, 19668, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 38418, 0, 3,
                                                                       36078, 17994, 36438, 7428,
                                                                       7563, 19938, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 38868, 0, 3,
                                                                       36438, 18210, 36798, 7563,
                                                                       7698, 20208, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 39318, 0, 3,
                                                                       36798, 18426, 37158, 7698,
                                                                       7833, 20478, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 39768, 0, 3,
                                                                       37518, 19398, 37968, 8103,
                                                                       8268, 21408, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 40318, 0, 3,
                                                                       37968, 19668, 38418, 8268,
                                                                       8433, 21738, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 40868, 0, 3,
                                                                       38418, 19938, 38868, 8433,
                                                                       8598, 22068, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 41418, 0, 3,
                                                                       38868, 20208, 39318, 8598,
                                                                       8763, 22398, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 41968, 0, 3,
                                                                       39768, 21408, 40318, 9093,
                                                                       9291, 23520, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 42628, 0, 3,
                                                                       40318, 21738, 40868, 9291,
                                                                       9489, 23916, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 43288, 0, 3,
                                                                       40868, 22068, 41418, 9489,
                                                                       9687, 24312, ncols, gamma,
                                                                       p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 43948, 0, 3,
                                                                       41968, 23520, 42628,
                                                                       10083, 10317, 25644,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 44728, 0, 3,
                                                                       42628, 23916, 43288,
                                                                       10317, 10551, 26112,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 45508, 0, 3,
                                                                       43948, 25644, 44728,
                                                                       11019, 11292, 27672,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46418, 3, 11838,
                                                                       11844, 28218, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46433, 3, 11844,
                                                                       11850, 28228, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46448, 3, 11850,
                                                                       11856, 28238, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46463, 3, 11856,
                                                                       11862, 28248, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46478, 3, 11862,
                                                                       11868, 28258, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46493, 3, 11868,
                                                                       11874, 28268, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46508, 3, 11874,
                                                                       11880, 28278, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46523, 3, 11880,
                                                                       11886, 28288, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46538, 3, 11886,
                                                                       11892, 28298, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46553, 3, 11892,
                                                                       11898, 28308, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46568, 3, 11898,
                                                                       11904, 28318, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46583, 3, 11904,
                                                                       11910, 28328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46598, 3, 11910,
                                                                       11916, 28338, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46613, 0, 3,
                                                                       46418, 28218, 46433,
                                                                       11928, 11946, 28348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46658, 0, 3,
                                                                       46433, 28228, 46448,
                                                                       11946, 11964, 28378,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46703, 0, 3,
                                                                       46448, 28238, 46463,
                                                                       11964, 11982, 28408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46748, 0, 3,
                                                                       46463, 28248, 46478,
                                                                       11982, 12000, 28438,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46793, 0, 3,
                                                                       46478, 28258, 46493,
                                                                       12000, 12018, 28468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46838, 0, 3,
                                                                       46493, 28268, 46508,
                                                                       12018, 12036, 28498,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46883, 0, 3,
                                                                       46508, 28278, 46523,
                                                                       12036, 12054, 28528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46928, 0, 3,
                                                                       46523, 28288, 46538,
                                                                       12054, 12072, 28558,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46973, 0, 3,
                                                                       46538, 28298, 46553,
                                                                       12072, 12090, 28588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 47018, 0, 3,
                                                                       46553, 28308, 46568,
                                                                       12090, 12108, 28618,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 47063, 0, 3,
                                                                       46568, 28318, 46583,
                                                                       12108, 12126, 28648,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 47108, 0, 3,
                                                                       46583, 28328, 46598,
                                                                       12126, 12144, 28678,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47153, 0, 3,
                                                                       46613, 28348, 46658,
                                                                       12180, 12216, 28708,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47243, 0, 3,
                                                                       46658, 28378, 46703,
                                                                       12216, 12252, 28768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47333, 0, 3,
                                                                       46703, 28408, 46748,
                                                                       12252, 12288, 28828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47423, 0, 3,
                                                                       46748, 28438, 46793,
                                                                       12288, 12324, 28888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47513, 0, 3,
                                                                       46793, 28468, 46838,
                                                                       12324, 12360, 28948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47603, 0, 3,
                                                                       46838, 28498, 46883,
                                                                       12360, 12396, 29008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47693, 0, 3,
                                                                       46883, 28528, 46928,
                                                                       12396, 12432, 29068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47783, 0, 3,
                                                                       46928, 28558, 46973,
                                                                       12432, 12468, 29128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47873, 0, 3,
                                                                       46973, 28588, 47018,
                                                                       12468, 12504, 29188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47963, 0, 3,
                                                                       47018, 28618, 47063,
                                                                       12504, 12540, 29248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 48053, 0, 3,
                                                                       47063, 28648, 47108,
                                                                       12540, 12576, 29308,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48143, 0, 3,
                                                                       47153, 28708, 47243,
                                                                       12648, 12708, 29368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48293, 0, 3,
                                                                       47243, 28768, 47333,
                                                                       12708, 12768, 29468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48443, 0, 3,
                                                                       47333, 28828, 47423,
                                                                       12768, 12828, 29568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48593, 0, 3,
                                                                       47423, 28888, 47513,
                                                                       12828, 12888, 29668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48743, 0, 3,
                                                                       47513, 28948, 47603,
                                                                       12888, 12948, 29768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48893, 0, 3,
                                                                       47603, 29008, 47693,
                                                                       12948, 13008, 29868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 49043, 0, 3,
                                                                       47693, 29068, 47783,
                                                                       13008, 13068, 29968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 49193, 0, 3,
                                                                       47783, 29128, 47873,
                                                                       13068, 13128, 30068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 49343, 0, 3,
                                                                       47873, 29188, 47963,
                                                                       13128, 13188, 30168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 49493, 0, 3,
                                                                       47963, 29248, 48053,
                                                                       13188, 13248, 30268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49643, 0, 3,
                                                                       48143, 29368, 48293,
                                                                       13368, 13458, 30368,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49868, 0, 3,
                                                                       48293, 29468, 48443,
                                                                       13458, 13548, 30518,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50093, 0, 3,
                                                                       48443, 29568, 48593,
                                                                       13548, 13638, 30668,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50318, 0, 3,
                                                                       48593, 29668, 48743,
                                                                       13638, 13728, 30818,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50543, 0, 3,
                                                                       48743, 29768, 48893,
                                                                       13728, 13818, 30968,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50768, 0, 3,
                                                                       48893, 29868, 49043,
                                                                       13818, 13908, 31118,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50993, 0, 3,
                                                                       49043, 29968, 49193,
                                                                       13908, 13998, 31268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 51218, 0, 3,
                                                                       49193, 30068, 49343,
                                                                       13998, 14088, 31418,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 51443, 0, 3,
                                                                       49343, 30168, 49493,
                                                                       14088, 14178, 31568,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51668, 0, 3,
                                                                       49643, 30368, 49868,
                                                                       14358, 14484, 31718,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51983, 0, 3,
                                                                       49868, 30518, 50093,
                                                                       14484, 14610, 31928,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52298, 0, 3,
                                                                       50093, 30668, 50318,
                                                                       14610, 14736, 32138,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52613, 0, 3,
                                                                       50318, 30818, 50543,
                                                                       14736, 14862, 32348,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52928, 0, 3,
                                                                       50543, 30968, 50768,
                                                                       14862, 14988, 32558,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 53243, 0, 3,
                                                                       50768, 31118, 50993,
                                                                       14988, 15114, 32768,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 53558, 0, 3,
                                                                       50993, 31268, 51218,
                                                                       15114, 15240, 32978,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 53873, 0, 3,
                                                                       51218, 31418, 51443,
                                                                       15240, 15366, 33188,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 54188, 0, 3,
                                                                       51668, 31718, 51983,
                                                                       15618, 15786, 33398,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 54608, 0, 3,
                                                                       51983, 31928, 52298,
                                                                       15786, 15954, 33678,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55028, 0, 3,
                                                                       52298, 32138, 52613,
                                                                       15954, 16122, 33958,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55448, 0, 3,
                                                                       52613, 32348, 52928,
                                                                       16122, 16290, 34238,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55868, 0, 3,
                                                                       52928, 32558, 53243,
                                                                       16290, 16458, 34518,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 56288, 0, 3,
                                                                       53243, 32768, 53558,
                                                                       16458, 16626, 34798,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 56708, 0, 3,
                                                                       53558, 32978, 53873,
                                                                       16626, 16794, 35078,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 57128, 0, 3,
                                                                       54188, 33398, 54608,
                                                                       17130, 17346, 35358,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 57668, 0, 3,
                                                                       54608, 33678, 55028,
                                                                       17346, 17562, 35718,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 58208, 0, 3,
                                                                       55028, 33958, 55448,
                                                                       17562, 17778, 36078,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 58748, 0, 3,
                                                                       55448, 34238, 55868,
                                                                       17778, 17994, 36438,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 59288, 0, 3,
                                                                       55868, 34518, 56288,
                                                                       17994, 18210, 36798,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 59828, 0, 3,
                                                                       56288, 34798, 56708,
                                                                       18210, 18426, 37158,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 60368, 0, 3,
                                                                       57128, 35358, 57668,
                                                                       18858, 19128, 37518,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 61043, 0, 3,
                                                                       57668, 35718, 58208,
                                                                       19128, 19398, 37968,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 61718, 0, 3,
                                                                       58208, 36078, 58748,
                                                                       19398, 19668, 38418,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 62393, 0, 3,
                                                                       58748, 36438, 59288,
                                                                       19668, 19938, 38868,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 63068, 0, 3,
                                                                       59288, 36798, 59828,
                                                                       19938, 20208, 39318,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 63743, 0, 3,
                                                                       60368, 37518, 61043,
                                                                       20748, 21078, 39768,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 64568, 0, 3,
                                                                       61043, 37968, 61718,
                                                                       21078, 21408, 40318,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 65393, 0, 3,
                                                                       61718, 38418, 62393,
                                                                       21408, 21738, 40868,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 66218, 0, 3,
                                                                       62393, 38868, 63068,
                                                                       21738, 22068, 41418,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 67043, 0, 3,
                                                                       63743, 39768, 64568,
                                                                       22728, 23124, 41968,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 68033, 0, 3,
                                                                       64568, 40318, 65393,
                                                                       23124, 23520, 42628,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 69023, 0, 3,
                                                                       65393, 40868, 66218,
                                                                       23520, 23916, 43288,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 70013, 0, 3,
                                                                       67043, 41968, 68033,
                                                                       24708, 25176, 43948,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 71183, 0, 3,
                                                                       68033, 42628, 69023,
                                                                       25176, 25644, 44728,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 72353, 0, 3,
                                                                       70013, 43948, 71183,
                                                                       26580, 27126, 45508,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 73718, 54188, 420, ncols);

                    simdfunc::contract_primitives(buffer, 74390, 57128, 540, ncols);

                    simdfunc::contract_primitives(buffer, 75254, 60368, 675, ncols);

                    simdfunc::contract_primitives(buffer, 76334, 63743, 825, ncols);

                    simdfunc::contract_primitives(buffer, 77654, 67043, 990, ncols);

                    simdfunc::contract_primitives(buffer, 79238, 70013, 1170, ncols);

                    simdfunc::contract_primitives(buffer, 81110, 72353, 1365, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 74138, 73718, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 74930, 74390, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 75929, 75254, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 77159, 76334, 55, 1, nmax);

        simdtrf::transform_g_inner(buffer, 78644, 77654, 66, 1, nmax);

        simdtrf::transform_g_inner(buffer, 80408, 79238, 78, 1, nmax);

        simdtrf::transform_g_inner(buffer, 82475, 81110, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 83294, 74138, 74930, 9, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 84050, 74930, 75929, 9, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 85022, 75929, 77159, 9, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 86237, 77159, 78644, 9, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 87722, 78644, 80408, 9, nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 89504, 80408, 82475, 9, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 91610, 83294, 84050, 9, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 93122, 84050, 85022, 9, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 95066, 85022, 86237, 9, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 97496, 86237, 87722, 9, nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 100466, 87722, 89504, 9,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 104030, 91610, 93122, 9,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 106550, 93122, 95066, 9,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 109790, 95066, 97496, 9,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 113840, 97496, 100466, 9,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 118790, 104030, 106550, 9,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 122570, 106550, 109790, 9,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 127430, 109790, 113840, 9,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 133505, 118790, 122570, 9,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 138797, 122570, 127430, 9,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 145601, 133505, 138797, 9,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 152657, 145601, 28, 9, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 152657, 117, nmax);
    }

    for (size_t m = 0; m < 1521; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
