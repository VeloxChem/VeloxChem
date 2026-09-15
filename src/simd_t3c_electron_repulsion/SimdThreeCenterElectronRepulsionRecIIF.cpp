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


#include "SimdThreeCenterElectronRepulsionRecIIF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
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
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iif_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iif_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 103328, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1183 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 103328, 40048, 6146, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
                                                        ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 8, 9,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 9, 10,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 11, 12,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 12, 13,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 13, 14,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 101, 0, 3, 14, 15,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 107, 0, 3, 15, 16,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 16, 17,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 119, 0, 3, 17, 18,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 125, 0, 3, 18, 19,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 131, 0, 3, 19, 20,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 137, 0, 3, 20, 21,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 23, 26,
                                                                       65, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 26, 29,
                                                                       71, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 29, 32,
                                                                       77, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 32, 35,
                                                                       83, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 35, 38,
                                                                       89, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 38, 41,
                                                                       95, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 41, 44,
                                                                       101, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 44, 47,
                                                                       107, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 47, 50,
                                                                       113, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 50, 53,
                                                                       119, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 53, 56,
                                                                       125, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 56, 59,
                                                                       131, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 263, 0, 3, 65, 71,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 71, 77,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 293, 0, 3, 77, 83,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 83, 89,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 323, 0, 3, 89, 95,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 95,
                                                                       101, 193, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 353, 0, 3, 101,
                                                                       107, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 107,
                                                                       113, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 383, 0, 3, 113,
                                                                       119, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 119,
                                                                       125, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 413, 0, 3, 125,
                                                                       131, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 143,
                                                                       153, 263, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 449, 0, 3, 153,
                                                                       163, 278, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 470, 0, 3, 163,
                                                                       173, 293, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 491, 0, 3, 173,
                                                                       183, 308, 323, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 512, 0, 3, 183,
                                                                       193, 323, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 533, 0, 3, 193,
                                                                       203, 338, 353, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 554, 0, 3, 203,
                                                                       213, 353, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 575, 0, 3, 213,
                                                                       223, 368, 383, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 596, 0, 3, 223,
                                                                       233, 383, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 617, 0, 3, 233,
                                                                       243, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 263,
                                                                       278, 428, 449, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 666, 0, 3, 278,
                                                                       293, 449, 470, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 694, 0, 3, 293,
                                                                       308, 470, 491, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 722, 0, 3, 308,
                                                                       323, 491, 512, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 750, 0, 3, 323,
                                                                       338, 512, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 778, 0, 3, 338,
                                                                       353, 533, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 806, 0, 3, 353,
                                                                       368, 554, 575, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 834, 0, 3, 368,
                                                                       383, 575, 596, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 862, 0, 3, 383,
                                                                       398, 596, 617, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 890, 0, 3, 428,
                                                                       449, 638, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 926, 0, 3, 449,
                                                                       470, 666, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 962, 0, 3, 470,
                                                                       491, 694, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 998, 0, 3, 491,
                                                                       512, 722, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1034, 0, 3, 512,
                                                                       533, 750, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1070, 0, 3, 533,
                                                                       554, 778, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1106, 0, 3, 554,
                                                                       575, 806, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 575,
                                                                       596, 834, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1178, 0, 3, 638,
                                                                       666, 890, 926, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1223, 0, 3, 666,
                                                                       694, 926, 962, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 694,
                                                                       722, 962, 998, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1313, 0, 3, 722,
                                                                       750, 998, 1034, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1358, 0, 3, 750,
                                                                       778, 1034, 1070, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1403, 0, 3, 778,
                                                                       806, 1070, 1106, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1448, 0, 3, 806,
                                                                       834, 1106, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1493, 0, 3, 890,
                                                                       926, 1178, 1223, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 926,
                                                                       962, 1223, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1603, 0, 3, 962,
                                                                       998, 1268, 1313, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1658, 0, 3, 998,
                                                                       1034, 1313, 1358, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1713, 0, 3, 1034,
                                                                       1070, 1358, 1403, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1768, 0, 3, 1070,
                                                                       1106, 1403, 1448, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1823, 0, 3, 1178,
                                                                       1223, 1493, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1889, 0, 3, 1223,
                                                                       1268, 1548, 1603, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1955, 0, 3, 1268,
                                                                       1313, 1603, 1658, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2021, 0, 3, 1313,
                                                                       1358, 1658, 1713, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2087, 0, 3, 1358,
                                                                       1403, 1713, 1768, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1493,
                                                                       1548, 1823, 1889, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2231, 0, 3, 1548,
                                                                       1603, 1889, 1955, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2309, 0, 3, 1603,
                                                                       1658, 1955, 2021, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2387, 0, 3, 1658,
                                                                       1713, 2021, 2087, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 2465, 0, 3, 1823,
                                                                       1889, 2153, 2231, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 2556, 0, 3, 1889,
                                                                       1955, 2231, 2309, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 2647, 0, 3, 1955,
                                                                       2021, 2309, 2387, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2738, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2741, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2744, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2747, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2750, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2753, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2756, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2759, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2762, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2765, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2768, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2771, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2774, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2777, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2780, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2783, 3, 8, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2792, 3, 9, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2801, 3, 10, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2810, 3, 11, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2819, 3, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2828, 3, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2837, 3, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2846, 3, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2855, 3, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2864, 3, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2873, 3, 18, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2882, 3, 19, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2891, 3, 20, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2900, 3, 21, 62,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2909, 3, 23, 65,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2927, 3, 26, 71,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2945, 3, 29, 77,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2963, 3, 32, 83,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2981, 3, 35, 89,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2999, 3, 38, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3017, 3, 41, 101,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3035, 3, 44, 107,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3053, 3, 47, 113,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3071, 3, 50, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3089, 3, 53, 125,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3107, 3, 56, 131,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3125, 3, 59, 137,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3143, 3, 65, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3173, 3, 71, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3203, 3, 77, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3233, 3, 83, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3263, 3, 89, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3293, 3, 95, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3323, 3, 101, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3353, 3, 107, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3383, 3, 113, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3413, 3, 119, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3443, 3, 125, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3473, 3, 131, 253,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3503, 3, 143, 263,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3548, 3, 153, 278,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3593, 3, 163, 293,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3638, 3, 173, 308,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3683, 3, 183, 323,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3728, 3, 193, 338,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3773, 3, 203, 353,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3818, 3, 213, 368,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3863, 3, 223, 383,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3908, 3, 233, 398,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3953, 3, 243, 413,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3998, 3, 263, 428,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4061, 3, 278, 449,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4124, 3, 293, 470,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4187, 3, 308, 491,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4250, 3, 323, 512,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4313, 3, 338, 533,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4376, 3, 353, 554,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4439, 3, 368, 575,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4502, 3, 383, 596,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4565, 3, 398, 617,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4628, 3, 428, 638,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4712, 3, 449, 666,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4796, 3, 470, 694,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4880, 3, 491, 722,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4964, 3, 512, 750,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5048, 3, 533, 778,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5132, 3, 554, 806,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5216, 3, 575, 834,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5300, 3, 596, 862,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5384, 3, 638, 890,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5492, 3, 666, 926,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5600, 3, 694, 962,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5708, 3, 722, 998,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5816, 3, 750,
                                                                       1034, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5924, 3, 778,
                                                                       1070, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6032, 3, 806,
                                                                       1106, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6140, 3, 834,
                                                                       1142, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6248, 3, 890,
                                                                       1178, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6383, 3, 926,
                                                                       1223, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6518, 3, 962,
                                                                       1268, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6653, 3, 998,
                                                                       1313, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6788, 3, 1034,
                                                                       1358, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6923, 3, 1070,
                                                                       1403, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7058, 3, 1106,
                                                                       1448, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7193, 3, 1178,
                                                                       1493, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7358, 3, 1223,
                                                                       1548, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7523, 3, 1268,
                                                                       1603, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7688, 3, 1313,
                                                                       1658, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7853, 3, 1358,
                                                                       1713, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8018, 3, 1403,
                                                                       1768, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8183, 3, 1493,
                                                                       1823, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8381, 3, 1548,
                                                                       1889, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8579, 3, 1603,
                                                                       1955, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8777, 3, 1658,
                                                                       2021, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8975, 3, 1713,
                                                                       2087, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 9173, 3, 1823,
                                                                       2153, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 9407, 3, 1889,
                                                                       2231, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 9641, 3, 1955,
                                                                       2309, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 9875, 3, 2021,
                                                                       2387, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 10109, 3, 2153,
                                                                       2465, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 10382, 3, 2231,
                                                                       2556, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 10655, 3, 2309,
                                                                       2647, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10928, 3, 8, 9,
                                                                       2744, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10934, 3, 9, 10,
                                                                       2747, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10940, 3, 10, 11,
                                                                       2750, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10946, 3, 11, 12,
                                                                       2753, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10952, 3, 12, 13,
                                                                       2756, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10958, 3, 13, 14,
                                                                       2759, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10964, 3, 14, 15,
                                                                       2762, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10970, 3, 15, 16,
                                                                       2765, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10976, 3, 16, 17,
                                                                       2768, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10982, 3, 17, 18,
                                                                       2771, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10988, 3, 18, 19,
                                                                       2774, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10994, 3, 19, 20,
                                                                       2777, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11000, 3, 20, 21,
                                                                       2780, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11006, 0, 3,
                                                                       10928, 2744, 10934, 2801,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11024, 0, 3,
                                                                       10934, 2747, 10940, 2810,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11042, 0, 3,
                                                                       10940, 2750, 10946, 2819,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11060, 0, 3,
                                                                       10946, 2753, 10952, 2828,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11078, 0, 3,
                                                                       10952, 2756, 10958, 2837,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11096, 0, 3,
                                                                       10958, 2759, 10964, 2846,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11114, 0, 3,
                                                                       10964, 2762, 10970, 2855,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11132, 0, 3,
                                                                       10970, 2765, 10976, 2864,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11150, 0, 3,
                                                                       10976, 2768, 10982, 2873,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11168, 0, 3,
                                                                       10982, 2771, 10988, 2882,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11186, 0, 3,
                                                                       10988, 2774, 10994, 2891,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11204, 0, 3,
                                                                       10994, 2777, 11000, 2900,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11222, 0, 3,
                                                                       11006, 2801, 11024, 65,
                                                                       71, 2945, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11258, 0, 3,
                                                                       11024, 2810, 11042, 71,
                                                                       77, 2963, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11294, 0, 3,
                                                                       11042, 2819, 11060, 77,
                                                                       83, 2981, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11330, 0, 3,
                                                                       11060, 2828, 11078, 83,
                                                                       89, 2999, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11366, 0, 3,
                                                                       11078, 2837, 11096, 89,
                                                                       95, 3017, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11402, 0, 3,
                                                                       11096, 2846, 11114, 95,
                                                                       101, 3035, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11438, 0, 3,
                                                                       11114, 2855, 11132, 101,
                                                                       107, 3053, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11474, 0, 3,
                                                                       11132, 2864, 11150, 107,
                                                                       113, 3071, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11510, 0, 3,
                                                                       11150, 2873, 11168, 113,
                                                                       119, 3089, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11546, 0, 3,
                                                                       11168, 2882, 11186, 119,
                                                                       125, 3107, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11582, 0, 3,
                                                                       11186, 2891, 11204, 125,
                                                                       131, 3125, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11618, 0, 3,
                                                                       11222, 2945, 11258, 143,
                                                                       153, 3203, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11678, 0, 3,
                                                                       11258, 2963, 11294, 153,
                                                                       163, 3233, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11738, 0, 3,
                                                                       11294, 2981, 11330, 163,
                                                                       173, 3263, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11798, 0, 3,
                                                                       11330, 2999, 11366, 173,
                                                                       183, 3293, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11858, 0, 3,
                                                                       11366, 3017, 11402, 183,
                                                                       193, 3323, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11918, 0, 3,
                                                                       11402, 3035, 11438, 193,
                                                                       203, 3353, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11978, 0, 3,
                                                                       11438, 3053, 11474, 203,
                                                                       213, 3383, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12038, 0, 3,
                                                                       11474, 3071, 11510, 213,
                                                                       223, 3413, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12098, 0, 3,
                                                                       11510, 3089, 11546, 223,
                                                                       233, 3443, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12158, 0, 3,
                                                                       11546, 3107, 11582, 233,
                                                                       243, 3473, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12218, 0, 3,
                                                                       11618, 3203, 11678, 263,
                                                                       278, 3593, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12308, 0, 3,
                                                                       11678, 3233, 11738, 278,
                                                                       293, 3638, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12398, 0, 3,
                                                                       11738, 3263, 11798, 293,
                                                                       308, 3683, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12488, 0, 3,
                                                                       11798, 3293, 11858, 308,
                                                                       323, 3728, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12578, 0, 3,
                                                                       11858, 3323, 11918, 323,
                                                                       338, 3773, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12668, 0, 3,
                                                                       11918, 3353, 11978, 338,
                                                                       353, 3818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12758, 0, 3,
                                                                       11978, 3383, 12038, 353,
                                                                       368, 3863, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12848, 0, 3,
                                                                       12038, 3413, 12098, 368,
                                                                       383, 3908, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12938, 0, 3,
                                                                       12098, 3443, 12158, 383,
                                                                       398, 3953, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13028, 0, 3,
                                                                       12218, 3593, 12308, 428,
                                                                       449, 4124, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13154, 0, 3,
                                                                       12308, 3638, 12398, 449,
                                                                       470, 4187, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13280, 0, 3,
                                                                       12398, 3683, 12488, 470,
                                                                       491, 4250, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13406, 0, 3,
                                                                       12488, 3728, 12578, 491,
                                                                       512, 4313, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13532, 0, 3,
                                                                       12578, 3773, 12668, 512,
                                                                       533, 4376, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13658, 0, 3,
                                                                       12668, 3818, 12758, 533,
                                                                       554, 4439, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13784, 0, 3,
                                                                       12758, 3863, 12848, 554,
                                                                       575, 4502, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13910, 0, 3,
                                                                       12848, 3908, 12938, 575,
                                                                       596, 4565, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14036, 0, 3,
                                                                       13028, 4124, 13154, 638,
                                                                       666, 4796, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14204, 0, 3,
                                                                       13154, 4187, 13280, 666,
                                                                       694, 4880, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14372, 0, 3,
                                                                       13280, 4250, 13406, 694,
                                                                       722, 4964, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14540, 0, 3,
                                                                       13406, 4313, 13532, 722,
                                                                       750, 5048, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14708, 0, 3,
                                                                       13532, 4376, 13658, 750,
                                                                       778, 5132, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14876, 0, 3,
                                                                       13658, 4439, 13784, 778,
                                                                       806, 5216, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15044, 0, 3,
                                                                       13784, 4502, 13910, 806,
                                                                       834, 5300, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15212, 0, 3,
                                                                       14036, 4796, 14204, 890,
                                                                       926, 5600, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15428, 0, 3,
                                                                       14204, 4880, 14372, 926,
                                                                       962, 5708, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15644, 0, 3,
                                                                       14372, 4964, 14540, 962,
                                                                       998, 5816, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15860, 0, 3,
                                                                       14540, 5048, 14708, 998,
                                                                       1034, 5924, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16076, 0, 3,
                                                                       14708, 5132, 14876, 1034,
                                                                       1070, 6032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16292, 0, 3,
                                                                       14876, 5216, 15044, 1070,
                                                                       1106, 6140, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16508, 0, 3,
                                                                       15212, 5600, 15428, 1178,
                                                                       1223, 6518, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16778, 0, 3,
                                                                       15428, 5708, 15644, 1223,
                                                                       1268, 6653, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17048, 0, 3,
                                                                       15644, 5816, 15860, 1268,
                                                                       1313, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17318, 0, 3,
                                                                       15860, 5924, 16076, 1313,
                                                                       1358, 6923, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17588, 0, 3,
                                                                       16076, 6032, 16292, 1358,
                                                                       1403, 7058, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 17858, 0, 3,
                                                                       16508, 6518, 16778, 1493,
                                                                       1548, 7523, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 18188, 0, 3,
                                                                       16778, 6653, 17048, 1548,
                                                                       1603, 7688, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 18518, 0, 3,
                                                                       17048, 6788, 17318, 1603,
                                                                       1658, 7853, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 18848, 0, 3,
                                                                       17318, 6923, 17588, 1658,
                                                                       1713, 8018, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 19178, 0, 3,
                                                                       17858, 7523, 18188, 1823,
                                                                       1889, 8579, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 19574, 0, 3,
                                                                       18188, 7688, 18518, 1889,
                                                                       1955, 8777, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 19970, 0, 3,
                                                                       18518, 7853, 18848, 1955,
                                                                       2021, 8975, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 20366, 0, 3,
                                                                       19178, 8579, 19574, 2153,
                                                                       2231, 9641, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 20834, 0, 3,
                                                                       19574, 8777, 19970, 2231,
                                                                       2309, 9875, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 21302, 0, 3,
                                                                       20366, 9641, 20834, 2465,
                                                                       2556, 10655, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21848, 3, 2738,
                                                                       2741, 10928, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21858, 3, 2741,
                                                                       2744, 10934, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21868, 3, 2744,
                                                                       2747, 10940, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21878, 3, 2747,
                                                                       2750, 10946, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21888, 3, 2750,
                                                                       2753, 10952, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21898, 3, 2753,
                                                                       2756, 10958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21908, 3, 2756,
                                                                       2759, 10964, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21918, 3, 2759,
                                                                       2762, 10970, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21928, 3, 2762,
                                                                       2765, 10976, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21938, 3, 2765,
                                                                       2768, 10982, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21948, 3, 2768,
                                                                       2771, 10988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21958, 3, 2771,
                                                                       2774, 10994, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21968, 3, 2774,
                                                                       2777, 11000, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21978, 0, 3,
                                                                       21848, 10928, 21858, 2783,
                                                                       2792, 11006, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22008, 0, 3,
                                                                       21858, 10934, 21868, 2792,
                                                                       2801, 11024, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22038, 0, 3,
                                                                       21868, 10940, 21878, 2801,
                                                                       2810, 11042, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22068, 0, 3,
                                                                       21878, 10946, 21888, 2810,
                                                                       2819, 11060, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22098, 0, 3,
                                                                       21888, 10952, 21898, 2819,
                                                                       2828, 11078, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22128, 0, 3,
                                                                       21898, 10958, 21908, 2828,
                                                                       2837, 11096, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22158, 0, 3,
                                                                       21908, 10964, 21918, 2837,
                                                                       2846, 11114, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22188, 0, 3,
                                                                       21918, 10970, 21928, 2846,
                                                                       2855, 11132, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22218, 0, 3,
                                                                       21928, 10976, 21938, 2855,
                                                                       2864, 11150, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22248, 0, 3,
                                                                       21938, 10982, 21948, 2864,
                                                                       2873, 11168, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22278, 0, 3,
                                                                       21948, 10988, 21958, 2873,
                                                                       2882, 11186, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22308, 0, 3,
                                                                       21958, 10994, 21968, 2882,
                                                                       2891, 11204, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22338, 0, 3,
                                                                       21978, 11006, 22008, 2909,
                                                                       2927, 11222, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22398, 0, 3,
                                                                       22008, 11024, 22038, 2927,
                                                                       2945, 11258, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22458, 0, 3,
                                                                       22038, 11042, 22068, 2945,
                                                                       2963, 11294, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22518, 0, 3,
                                                                       22068, 11060, 22098, 2963,
                                                                       2981, 11330, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22578, 0, 3,
                                                                       22098, 11078, 22128, 2981,
                                                                       2999, 11366, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22638, 0, 3,
                                                                       22128, 11096, 22158, 2999,
                                                                       3017, 11402, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22698, 0, 3,
                                                                       22158, 11114, 22188, 3017,
                                                                       3035, 11438, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22758, 0, 3,
                                                                       22188, 11132, 22218, 3035,
                                                                       3053, 11474, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22818, 0, 3,
                                                                       22218, 11150, 22248, 3053,
                                                                       3071, 11510, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22878, 0, 3,
                                                                       22248, 11168, 22278, 3071,
                                                                       3089, 11546, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22938, 0, 3,
                                                                       22278, 11186, 22308, 3089,
                                                                       3107, 11582, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22998, 0, 3,
                                                                       22338, 11222, 22398, 3143,
                                                                       3173, 11618, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23098, 0, 3,
                                                                       22398, 11258, 22458, 3173,
                                                                       3203, 11678, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23198, 0, 3,
                                                                       22458, 11294, 22518, 3203,
                                                                       3233, 11738, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23298, 0, 3,
                                                                       22518, 11330, 22578, 3233,
                                                                       3263, 11798, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23398, 0, 3,
                                                                       22578, 11366, 22638, 3263,
                                                                       3293, 11858, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23498, 0, 3,
                                                                       22638, 11402, 22698, 3293,
                                                                       3323, 11918, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23598, 0, 3,
                                                                       22698, 11438, 22758, 3323,
                                                                       3353, 11978, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23698, 0, 3,
                                                                       22758, 11474, 22818, 3353,
                                                                       3383, 12038, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23798, 0, 3,
                                                                       22818, 11510, 22878, 3383,
                                                                       3413, 12098, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23898, 0, 3,
                                                                       22878, 11546, 22938, 3413,
                                                                       3443, 12158, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23998, 0, 3,
                                                                       22998, 11618, 23098, 3503,
                                                                       3548, 12218, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24148, 0, 3,
                                                                       23098, 11678, 23198, 3548,
                                                                       3593, 12308, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24298, 0, 3,
                                                                       23198, 11738, 23298, 3593,
                                                                       3638, 12398, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24448, 0, 3,
                                                                       23298, 11798, 23398, 3638,
                                                                       3683, 12488, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24598, 0, 3,
                                                                       23398, 11858, 23498, 3683,
                                                                       3728, 12578, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24748, 0, 3,
                                                                       23498, 11918, 23598, 3728,
                                                                       3773, 12668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24898, 0, 3,
                                                                       23598, 11978, 23698, 3773,
                                                                       3818, 12758, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 25048, 0, 3,
                                                                       23698, 12038, 23798, 3818,
                                                                       3863, 12848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 25198, 0, 3,
                                                                       23798, 12098, 23898, 3863,
                                                                       3908, 12938, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25348, 0, 3,
                                                                       23998, 12218, 24148, 3998,
                                                                       4061, 13028, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25558, 0, 3,
                                                                       24148, 12308, 24298, 4061,
                                                                       4124, 13154, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25768, 0, 3,
                                                                       24298, 12398, 24448, 4124,
                                                                       4187, 13280, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25978, 0, 3,
                                                                       24448, 12488, 24598, 4187,
                                                                       4250, 13406, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26188, 0, 3,
                                                                       24598, 12578, 24748, 4250,
                                                                       4313, 13532, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26398, 0, 3,
                                                                       24748, 12668, 24898, 4313,
                                                                       4376, 13658, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26608, 0, 3,
                                                                       24898, 12758, 25048, 4376,
                                                                       4439, 13784, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26818, 0, 3,
                                                                       25048, 12848, 25198, 4439,
                                                                       4502, 13910, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27028, 0, 3,
                                                                       25348, 13028, 25558, 4628,
                                                                       4712, 14036, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27308, 0, 3,
                                                                       25558, 13154, 25768, 4712,
                                                                       4796, 14204, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27588, 0, 3,
                                                                       25768, 13280, 25978, 4796,
                                                                       4880, 14372, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27868, 0, 3,
                                                                       25978, 13406, 26188, 4880,
                                                                       4964, 14540, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28148, 0, 3,
                                                                       26188, 13532, 26398, 4964,
                                                                       5048, 14708, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28428, 0, 3,
                                                                       26398, 13658, 26608, 5048,
                                                                       5132, 14876, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28708, 0, 3,
                                                                       26608, 13784, 26818, 5132,
                                                                       5216, 15044, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 28988, 0, 3,
                                                                       27028, 14036, 27308, 5384,
                                                                       5492, 15212, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 29348, 0, 3,
                                                                       27308, 14204, 27588, 5492,
                                                                       5600, 15428, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 29708, 0, 3,
                                                                       27588, 14372, 27868, 5600,
                                                                       5708, 15644, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30068, 0, 3,
                                                                       27868, 14540, 28148, 5708,
                                                                       5816, 15860, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30428, 0, 3,
                                                                       28148, 14708, 28428, 5816,
                                                                       5924, 16076, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30788, 0, 3,
                                                                       28428, 14876, 28708, 5924,
                                                                       6032, 16292, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 31148, 0, 3,
                                                                       28988, 15212, 29348, 6248,
                                                                       6383, 16508, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 31598, 0, 3,
                                                                       29348, 15428, 29708, 6383,
                                                                       6518, 16778, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 32048, 0, 3,
                                                                       29708, 15644, 30068, 6518,
                                                                       6653, 17048, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 32498, 0, 3,
                                                                       30068, 15860, 30428, 6653,
                                                                       6788, 17318, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 32948, 0, 3,
                                                                       30428, 16076, 30788, 6788,
                                                                       6923, 17588, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 33398, 0, 3,
                                                                       31148, 16508, 31598, 7193,
                                                                       7358, 17858, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 33948, 0, 3,
                                                                       31598, 16778, 32048, 7358,
                                                                       7523, 18188, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 34498, 0, 3,
                                                                       32048, 17048, 32498, 7523,
                                                                       7688, 18518, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 35048, 0, 3,
                                                                       32498, 17318, 32948, 7688,
                                                                       7853, 18848, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 35598, 0, 3,
                                                                       33398, 17858, 33948, 8183,
                                                                       8381, 19178, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 36258, 0, 3,
                                                                       33948, 18188, 34498, 8381,
                                                                       8579, 19574, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 36918, 0, 3,
                                                                       34498, 18518, 35048, 8579,
                                                                       8777, 19970, ncols, gamma,
                                                                       p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 37578, 0, 3,
                                                                       35598, 19178, 36258, 9173,
                                                                       9407, 20366, ncols, gamma,
                                                                       p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 38358, 0, 3,
                                                                       36258, 19574, 36918, 9407,
                                                                       9641, 20834, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 39138, 0, 3,
                                                                       37578, 20366, 38358,
                                                                       10109, 10382, 21302,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 40048, 27028, 280, ncols);

                    simdfunc::contract_primitives(buffer, 40524, 28988, 360, ncols);

                    simdfunc::contract_primitives(buffer, 41136, 31148, 450, ncols);

                    simdfunc::contract_primitives(buffer, 41901, 33398, 550, ncols);

                    simdfunc::contract_primitives(buffer, 42836, 35598, 660, ncols);

                    simdfunc::contract_primitives(buffer, 43958, 37578, 780, ncols);

                    simdfunc::contract_primitives(buffer, 45284, 39138, 910, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 40328, 40048, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 40884, 40524, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 41586, 41136, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 42451, 41901, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 43496, 42836, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 44738, 43958, 78, 1, nmax);

        simdtrf::transform_f_inner(buffer, 46194, 45284, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 46831, 40328, 40884, 7, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 47419, 40884, 41586, 7, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 48175, 41586, 42451, 7, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 49120, 42451, 43496, 7, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 50275, 43496, 44738, 7, nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 51661, 44738, 46194, 7, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 53299, 46831, 47419, 7, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 54475, 47419, 48175, 7, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 55987, 48175, 49120, 7, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 57877, 49120, 50275, 7, nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 60187, 50275, 51661, 7, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 62959, 53299, 54475, 7, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 64919, 54475, 55987, 7, nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 67439, 55987, 57877, 7, nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 70589, 57877, 60187, 7, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 74439, 62959, 64919, 7, nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 77379, 64919, 67439, 7, nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 81159, 67439, 70589, 7, nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 85884, 74439, 77379, 7, nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 90000, 77379, 81159, 7, nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 95292, 85884, 90000, 7, nmax);

        simdtrf::transform_i_inner(buffer, 100780, 95292, 28, 7, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 100780, 91, nmax);
    }

    for (size_t m = 0; m < 1183; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
