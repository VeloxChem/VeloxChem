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


#include "SimdThreeCenterElectronRepulsionRecIHF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
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
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ihf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ihf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 67796, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1001 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 67796, 30220, 4690, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
                                                        ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 8, 9,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 9, 10,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 10, 11,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 11, 12,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 12, 13,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 91, 0, 3, 13, 14,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 14, 15,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 15, 16,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 109, 0, 3, 16, 17,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 115, 0, 3, 17, 18,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 121, 0, 3, 18, 19,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 127, 0, 3, 19, 20,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 22, 25,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 25, 28,
                                                                       67, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 28, 31,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 31, 34,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 34, 37,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 37, 40,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 40, 43,
                                                                       97, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 43, 46,
                                                                       103, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 46, 49,
                                                                       109, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 49, 52,
                                                                       115, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 52, 55,
                                                                       121, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 61, 67,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 67, 73,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 73, 79,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 79, 85,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 303, 0, 3, 85, 91,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 91, 97,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 333, 0, 3, 97,
                                                                       103, 193, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 103,
                                                                       109, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 363, 0, 3, 109,
                                                                       115, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 115,
                                                                       121, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 393, 0, 3, 133,
                                                                       143, 243, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 414, 0, 3, 143,
                                                                       153, 258, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 435, 0, 3, 153,
                                                                       163, 273, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 456, 0, 3, 163,
                                                                       173, 288, 303, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 477, 0, 3, 173,
                                                                       183, 303, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 183,
                                                                       193, 318, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 519, 0, 3, 193,
                                                                       203, 333, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 540, 0, 3, 203,
                                                                       213, 348, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 561, 0, 3, 213,
                                                                       223, 363, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 582, 0, 3, 243,
                                                                       258, 393, 414, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 610, 0, 3, 258,
                                                                       273, 414, 435, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 273,
                                                                       288, 435, 456, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 666, 0, 3, 288,
                                                                       303, 456, 477, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 694, 0, 3, 303,
                                                                       318, 477, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 722, 0, 3, 318,
                                                                       333, 498, 519, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 750, 0, 3, 333,
                                                                       348, 519, 540, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 778, 0, 3, 348,
                                                                       363, 540, 561, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 806, 0, 3, 393,
                                                                       414, 582, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 842, 0, 3, 414,
                                                                       435, 610, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 878, 0, 3, 435,
                                                                       456, 638, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 914, 0, 3, 456,
                                                                       477, 666, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 950, 0, 3, 477,
                                                                       498, 694, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 986, 0, 3, 498,
                                                                       519, 722, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1022, 0, 3, 519,
                                                                       540, 750, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 582,
                                                                       610, 806, 842, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1103, 0, 3, 610,
                                                                       638, 842, 878, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1148, 0, 3, 638,
                                                                       666, 878, 914, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1193, 0, 3, 666,
                                                                       694, 914, 950, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1238, 0, 3, 694,
                                                                       722, 950, 986, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1283, 0, 3, 722,
                                                                       750, 986, 1022, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1328, 0, 3, 806,
                                                                       842, 1058, 1103, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1383, 0, 3, 842,
                                                                       878, 1103, 1148, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1438, 0, 3, 878,
                                                                       914, 1148, 1193, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1493, 0, 3, 914,
                                                                       950, 1193, 1238, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 950,
                                                                       986, 1238, 1283, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1603, 0, 3, 1058,
                                                                       1103, 1328, 1383, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1669, 0, 3, 1103,
                                                                       1148, 1383, 1438, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1735, 0, 3, 1148,
                                                                       1193, 1438, 1493, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1801, 0, 3, 1193,
                                                                       1238, 1493, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 1867, 0, 3, 1328,
                                                                       1383, 1603, 1669, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 1945, 0, 3, 1383,
                                                                       1438, 1669, 1735, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2023, 0, 3, 1438,
                                                                       1493, 1735, 1801, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2101, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2104, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2107, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2110, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2113, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2116, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2119, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2122, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2125, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2128, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2131, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2134, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2137, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2140, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2143, 3, 8, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2152, 3, 9, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2161, 3, 10, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2170, 3, 11, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2179, 3, 12, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2188, 3, 13, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2197, 3, 14, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2206, 3, 15, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2215, 3, 16, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2224, 3, 17, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2233, 3, 18, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2242, 3, 19, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2251, 3, 20, 58,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2260, 3, 22, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2278, 3, 25, 67,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2296, 3, 28, 73,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2314, 3, 31, 79,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2332, 3, 34, 85,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2350, 3, 37, 91,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2368, 3, 40, 97,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2386, 3, 43, 103,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2404, 3, 46, 109,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2422, 3, 49, 115,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2440, 3, 52, 121,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2458, 3, 55, 127,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2476, 3, 61, 133,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2506, 3, 67, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2536, 3, 73, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2566, 3, 79, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2596, 3, 85, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2626, 3, 91, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2656, 3, 97, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2686, 3, 103, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2716, 3, 109, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2746, 3, 115, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2776, 3, 121, 233,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2806, 3, 133, 243,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2851, 3, 143, 258,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2896, 3, 153, 273,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2941, 3, 163, 288,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2986, 3, 173, 303,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3031, 3, 183, 318,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3076, 3, 193, 333,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3121, 3, 203, 348,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3166, 3, 213, 363,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3211, 3, 223, 378,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3256, 3, 243, 393,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3319, 3, 258, 414,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3382, 3, 273, 435,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3445, 3, 288, 456,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3508, 3, 303, 477,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3571, 3, 318, 498,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3634, 3, 333, 519,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3697, 3, 348, 540,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3760, 3, 363, 561,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3823, 3, 393, 582,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3907, 3, 414, 610,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3991, 3, 435, 638,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4075, 3, 456, 666,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4159, 3, 477, 694,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4243, 3, 498, 722,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4327, 3, 519, 750,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4411, 3, 540, 778,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4495, 3, 582, 806,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4603, 3, 610, 842,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4711, 3, 638, 878,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4819, 3, 666, 914,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4927, 3, 694, 950,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5035, 3, 722, 986,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5143, 3, 750,
                                                                       1022, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5251, 3, 806,
                                                                       1058, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5386, 3, 842,
                                                                       1103, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5521, 3, 878,
                                                                       1148, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5656, 3, 914,
                                                                       1193, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5791, 3, 950,
                                                                       1238, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5926, 3, 986,
                                                                       1283, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6061, 3, 1058,
                                                                       1328, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6226, 3, 1103,
                                                                       1383, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6391, 3, 1148,
                                                                       1438, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6556, 3, 1193,
                                                                       1493, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6721, 3, 1238,
                                                                       1548, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 6886, 3, 1328,
                                                                       1603, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 7084, 3, 1383,
                                                                       1669, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 7282, 3, 1438,
                                                                       1735, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 7480, 3, 1493,
                                                                       1801, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 7678, 3, 1603,
                                                                       1867, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 7912, 3, 1669,
                                                                       1945, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 8146, 3, 1735,
                                                                       2023, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8380, 3, 8, 9,
                                                                       2107, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8386, 3, 9, 10,
                                                                       2110, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8392, 3, 10, 11,
                                                                       2113, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8398, 3, 11, 12,
                                                                       2116, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8404, 3, 12, 13,
                                                                       2119, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8410, 3, 13, 14,
                                                                       2122, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8416, 3, 14, 15,
                                                                       2125, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8422, 3, 15, 16,
                                                                       2128, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8428, 3, 16, 17,
                                                                       2131, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8434, 3, 17, 18,
                                                                       2134, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8440, 3, 18, 19,
                                                                       2137, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8446, 3, 19, 20,
                                                                       2140, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8452, 0, 3, 8380,
                                                                       2107, 8386, 2161, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8470, 0, 3, 8386,
                                                                       2110, 8392, 2170, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8488, 0, 3, 8392,
                                                                       2113, 8398, 2179, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8506, 0, 3, 8398,
                                                                       2116, 8404, 2188, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8524, 0, 3, 8404,
                                                                       2119, 8410, 2197, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8542, 0, 3, 8410,
                                                                       2122, 8416, 2206, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8560, 0, 3, 8416,
                                                                       2125, 8422, 2215, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8578, 0, 3, 8422,
                                                                       2128, 8428, 2224, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8596, 0, 3, 8428,
                                                                       2131, 8434, 2233, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8614, 0, 3, 8434,
                                                                       2134, 8440, 2242, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8632, 0, 3, 8440,
                                                                       2137, 8446, 2251, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8650, 0, 3, 8452,
                                                                       2161, 8470, 61, 67, 2296,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8686, 0, 3, 8470,
                                                                       2170, 8488, 67, 73, 2314,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8722, 0, 3, 8488,
                                                                       2179, 8506, 73, 79, 2332,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8758, 0, 3, 8506,
                                                                       2188, 8524, 79, 85, 2350,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8794, 0, 3, 8524,
                                                                       2197, 8542, 85, 91, 2368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8830, 0, 3, 8542,
                                                                       2206, 8560, 91, 97, 2386,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8866, 0, 3, 8560,
                                                                       2215, 8578, 97, 103, 2404,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8902, 0, 3, 8578,
                                                                       2224, 8596, 103, 109,
                                                                       2422, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8938, 0, 3, 8596,
                                                                       2233, 8614, 109, 115,
                                                                       2440, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8974, 0, 3, 8614,
                                                                       2242, 8632, 115, 121,
                                                                       2458, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9010, 0, 3, 8650,
                                                                       2296, 8686, 133, 143,
                                                                       2536, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9070, 0, 3, 8686,
                                                                       2314, 8722, 143, 153,
                                                                       2566, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9130, 0, 3, 8722,
                                                                       2332, 8758, 153, 163,
                                                                       2596, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9190, 0, 3, 8758,
                                                                       2350, 8794, 163, 173,
                                                                       2626, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9250, 0, 3, 8794,
                                                                       2368, 8830, 173, 183,
                                                                       2656, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9310, 0, 3, 8830,
                                                                       2386, 8866, 183, 193,
                                                                       2686, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9370, 0, 3, 8866,
                                                                       2404, 8902, 193, 203,
                                                                       2716, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9430, 0, 3, 8902,
                                                                       2422, 8938, 203, 213,
                                                                       2746, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9490, 0, 3, 8938,
                                                                       2440, 8974, 213, 223,
                                                                       2776, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9550, 0, 3, 9010,
                                                                       2536, 9070, 243, 258,
                                                                       2896, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9640, 0, 3, 9070,
                                                                       2566, 9130, 258, 273,
                                                                       2941, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9730, 0, 3, 9130,
                                                                       2596, 9190, 273, 288,
                                                                       2986, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9820, 0, 3, 9190,
                                                                       2626, 9250, 288, 303,
                                                                       3031, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9910, 0, 3, 9250,
                                                                       2656, 9310, 303, 318,
                                                                       3076, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10000, 0, 3, 9310,
                                                                       2686, 9370, 318, 333,
                                                                       3121, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10090, 0, 3, 9370,
                                                                       2716, 9430, 333, 348,
                                                                       3166, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10180, 0, 3, 9430,
                                                                       2746, 9490, 348, 363,
                                                                       3211, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10270, 0, 3, 9550,
                                                                       2896, 9640, 393, 414,
                                                                       3382, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10396, 0, 3, 9640,
                                                                       2941, 9730, 414, 435,
                                                                       3445, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10522, 0, 3, 9730,
                                                                       2986, 9820, 435, 456,
                                                                       3508, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10648, 0, 3, 9820,
                                                                       3031, 9910, 456, 477,
                                                                       3571, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10774, 0, 3, 9910,
                                                                       3076, 10000, 477, 498,
                                                                       3634, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10900, 0, 3,
                                                                       10000, 3121, 10090, 498,
                                                                       519, 3697, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11026, 0, 3,
                                                                       10090, 3166, 10180, 519,
                                                                       540, 3760, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11152, 0, 3,
                                                                       10270, 3382, 10396, 582,
                                                                       610, 3991, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11320, 0, 3,
                                                                       10396, 3445, 10522, 610,
                                                                       638, 4075, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11488, 0, 3,
                                                                       10522, 3508, 10648, 638,
                                                                       666, 4159, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11656, 0, 3,
                                                                       10648, 3571, 10774, 666,
                                                                       694, 4243, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11824, 0, 3,
                                                                       10774, 3634, 10900, 694,
                                                                       722, 4327, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11992, 0, 3,
                                                                       10900, 3697, 11026, 722,
                                                                       750, 4411, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12160, 0, 3,
                                                                       11152, 3991, 11320, 806,
                                                                       842, 4711, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12376, 0, 3,
                                                                       11320, 4075, 11488, 842,
                                                                       878, 4819, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12592, 0, 3,
                                                                       11488, 4159, 11656, 878,
                                                                       914, 4927, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12808, 0, 3,
                                                                       11656, 4243, 11824, 914,
                                                                       950, 5035, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13024, 0, 3,
                                                                       11824, 4327, 11992, 950,
                                                                       986, 5143, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13240, 0, 3,
                                                                       12160, 4711, 12376, 1058,
                                                                       1103, 5521, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13510, 0, 3,
                                                                       12376, 4819, 12592, 1103,
                                                                       1148, 5656, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13780, 0, 3,
                                                                       12592, 4927, 12808, 1148,
                                                                       1193, 5791, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14050, 0, 3,
                                                                       12808, 5035, 13024, 1193,
                                                                       1238, 5926, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 14320, 0, 3,
                                                                       13240, 5521, 13510, 1328,
                                                                       1383, 6391, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 14650, 0, 3,
                                                                       13510, 5656, 13780, 1383,
                                                                       1438, 6556, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 14980, 0, 3,
                                                                       13780, 5791, 14050, 1438,
                                                                       1493, 6721, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 15310, 0, 3,
                                                                       14320, 6391, 14650, 1603,
                                                                       1669, 7282, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 15706, 0, 3,
                                                                       14650, 6556, 14980, 1669,
                                                                       1735, 7480, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 16102, 0, 3,
                                                                       15310, 7282, 15706, 1867,
                                                                       1945, 8146, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16570, 3, 2101,
                                                                       2104, 8380, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16580, 3, 2104,
                                                                       2107, 8386, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16590, 3, 2107,
                                                                       2110, 8392, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16600, 3, 2110,
                                                                       2113, 8398, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16610, 3, 2113,
                                                                       2116, 8404, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16620, 3, 2116,
                                                                       2119, 8410, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16630, 3, 2119,
                                                                       2122, 8416, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16640, 3, 2122,
                                                                       2125, 8422, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16650, 3, 2125,
                                                                       2128, 8428, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16660, 3, 2128,
                                                                       2131, 8434, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16670, 3, 2131,
                                                                       2134, 8440, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16680, 3, 2134,
                                                                       2137, 8446, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16690, 0, 3,
                                                                       16570, 8380, 16580, 2143,
                                                                       2152, 8452, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16720, 0, 3,
                                                                       16580, 8386, 16590, 2152,
                                                                       2161, 8470, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16750, 0, 3,
                                                                       16590, 8392, 16600, 2161,
                                                                       2170, 8488, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16780, 0, 3,
                                                                       16600, 8398, 16610, 2170,
                                                                       2179, 8506, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16810, 0, 3,
                                                                       16610, 8404, 16620, 2179,
                                                                       2188, 8524, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16840, 0, 3,
                                                                       16620, 8410, 16630, 2188,
                                                                       2197, 8542, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16870, 0, 3,
                                                                       16630, 8416, 16640, 2197,
                                                                       2206, 8560, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16900, 0, 3,
                                                                       16640, 8422, 16650, 2206,
                                                                       2215, 8578, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16930, 0, 3,
                                                                       16650, 8428, 16660, 2215,
                                                                       2224, 8596, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16960, 0, 3,
                                                                       16660, 8434, 16670, 2224,
                                                                       2233, 8614, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16990, 0, 3,
                                                                       16670, 8440, 16680, 2233,
                                                                       2242, 8632, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17020, 0, 3,
                                                                       16690, 8452, 16720, 2260,
                                                                       2278, 8650, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17080, 0, 3,
                                                                       16720, 8470, 16750, 2278,
                                                                       2296, 8686, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17140, 0, 3,
                                                                       16750, 8488, 16780, 2296,
                                                                       2314, 8722, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17200, 0, 3,
                                                                       16780, 8506, 16810, 2314,
                                                                       2332, 8758, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17260, 0, 3,
                                                                       16810, 8524, 16840, 2332,
                                                                       2350, 8794, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17320, 0, 3,
                                                                       16840, 8542, 16870, 2350,
                                                                       2368, 8830, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17380, 0, 3,
                                                                       16870, 8560, 16900, 2368,
                                                                       2386, 8866, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17440, 0, 3,
                                                                       16900, 8578, 16930, 2386,
                                                                       2404, 8902, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17500, 0, 3,
                                                                       16930, 8596, 16960, 2404,
                                                                       2422, 8938, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17560, 0, 3,
                                                                       16960, 8614, 16990, 2422,
                                                                       2440, 8974, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17620, 0, 3,
                                                                       17020, 8650, 17080, 2476,
                                                                       2506, 9010, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17720, 0, 3,
                                                                       17080, 8686, 17140, 2506,
                                                                       2536, 9070, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17820, 0, 3,
                                                                       17140, 8722, 17200, 2536,
                                                                       2566, 9130, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17920, 0, 3,
                                                                       17200, 8758, 17260, 2566,
                                                                       2596, 9190, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18020, 0, 3,
                                                                       17260, 8794, 17320, 2596,
                                                                       2626, 9250, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18120, 0, 3,
                                                                       17320, 8830, 17380, 2626,
                                                                       2656, 9310, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18220, 0, 3,
                                                                       17380, 8866, 17440, 2656,
                                                                       2686, 9370, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18320, 0, 3,
                                                                       17440, 8902, 17500, 2686,
                                                                       2716, 9430, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18420, 0, 3,
                                                                       17500, 8938, 17560, 2716,
                                                                       2746, 9490, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18520, 0, 3,
                                                                       17620, 9010, 17720, 2806,
                                                                       2851, 9550, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18670, 0, 3,
                                                                       17720, 9070, 17820, 2851,
                                                                       2896, 9640, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18820, 0, 3,
                                                                       17820, 9130, 17920, 2896,
                                                                       2941, 9730, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18970, 0, 3,
                                                                       17920, 9190, 18020, 2941,
                                                                       2986, 9820, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 19120, 0, 3,
                                                                       18020, 9250, 18120, 2986,
                                                                       3031, 9910, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 19270, 0, 3,
                                                                       18120, 9310, 18220, 3031,
                                                                       3076, 10000, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 19420, 0, 3,
                                                                       18220, 9370, 18320, 3076,
                                                                       3121, 10090, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 19570, 0, 3,
                                                                       18320, 9430, 18420, 3121,
                                                                       3166, 10180, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19720, 0, 3,
                                                                       18520, 9550, 18670, 3256,
                                                                       3319, 10270, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19930, 0, 3,
                                                                       18670, 9640, 18820, 3319,
                                                                       3382, 10396, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20140, 0, 3,
                                                                       18820, 9730, 18970, 3382,
                                                                       3445, 10522, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20350, 0, 3,
                                                                       18970, 9820, 19120, 3445,
                                                                       3508, 10648, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20560, 0, 3,
                                                                       19120, 9910, 19270, 3508,
                                                                       3571, 10774, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20770, 0, 3,
                                                                       19270, 10000, 19420, 3571,
                                                                       3634, 10900, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20980, 0, 3,
                                                                       19420, 10090, 19570, 3634,
                                                                       3697, 11026, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21190, 0, 3,
                                                                       19720, 10270, 19930, 3823,
                                                                       3907, 11152, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21470, 0, 3,
                                                                       19930, 10396, 20140, 3907,
                                                                       3991, 11320, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21750, 0, 3,
                                                                       20140, 10522, 20350, 3991,
                                                                       4075, 11488, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22030, 0, 3,
                                                                       20350, 10648, 20560, 4075,
                                                                       4159, 11656, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22310, 0, 3,
                                                                       20560, 10774, 20770, 4159,
                                                                       4243, 11824, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22590, 0, 3,
                                                                       20770, 10900, 20980, 4243,
                                                                       4327, 11992, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 22870, 0, 3,
                                                                       21190, 11152, 21470, 4495,
                                                                       4603, 12160, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23230, 0, 3,
                                                                       21470, 11320, 21750, 4603,
                                                                       4711, 12376, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23590, 0, 3,
                                                                       21750, 11488, 22030, 4711,
                                                                       4819, 12592, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23950, 0, 3,
                                                                       22030, 11656, 22310, 4819,
                                                                       4927, 12808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 24310, 0, 3,
                                                                       22310, 11824, 22590, 4927,
                                                                       5035, 13024, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 24670, 0, 3,
                                                                       22870, 12160, 23230, 5251,
                                                                       5386, 13240, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 25120, 0, 3,
                                                                       23230, 12376, 23590, 5386,
                                                                       5521, 13510, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 25570, 0, 3,
                                                                       23590, 12592, 23950, 5521,
                                                                       5656, 13780, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 26020, 0, 3,
                                                                       23950, 12808, 24310, 5656,
                                                                       5791, 14050, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 26470, 0, 3,
                                                                       24670, 13240, 25120, 6061,
                                                                       6226, 14320, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 27020, 0, 3,
                                                                       25120, 13510, 25570, 6226,
                                                                       6391, 14650, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 27570, 0, 3,
                                                                       25570, 13780, 26020, 6391,
                                                                       6556, 14980, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 28120, 0, 3,
                                                                       26470, 14320, 27020, 6886,
                                                                       7084, 15310, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 28780, 0, 3,
                                                                       27020, 14650, 27570, 7084,
                                                                       7282, 15706, ncols, gamma,
                                                                       p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 29440, 0, 3,
                                                                       28120, 15310, 28780, 7678,
                                                                       7912, 16102, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 30220, 21190, 280, ncols);

                    simdfunc::contract_primitives(buffer, 30696, 22870, 360, ncols);

                    simdfunc::contract_primitives(buffer, 31308, 24670, 450, ncols);

                    simdfunc::contract_primitives(buffer, 32073, 26470, 550, ncols);

                    simdfunc::contract_primitives(buffer, 33008, 28120, 660, ncols);

                    simdfunc::contract_primitives(buffer, 34130, 29440, 780, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 30500, 30220, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31056, 30696, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31758, 31308, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 32623, 32073, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 33668, 33008, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 34910, 34130, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 35456, 30500, 31056, 7, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 36044, 31056, 31758, 7, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 36800, 31758, 32623, 7, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 37745, 32623, 33668, 7, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 38900, 33668, 34910, 7, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 40286, 35456, 36044, 7, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 41462, 36044, 36800, 7, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 42974, 36800, 37745, 7, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 44864, 37745, 38900, 7, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 47174, 40286, 41462, 7, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 49134, 41462, 42974, 7, nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 51654, 42974, 44864, 7, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 54804, 47174, 49134, 7, nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 57744, 49134, 51654, 7, nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 61524, 54804, 57744, 7, nmax);

        simdtrf::transform_h_inner(buffer, 65640, 61524, 28, 7, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 65640, 77, nmax);
    }

    for (size_t m = 0; m < 1001; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
