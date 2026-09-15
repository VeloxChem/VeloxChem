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


#include "SimdThreeCenterElectronRepulsionRecGGK.hpp"

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
#include "SimdTransferGD.hpp"
#include "SimdTransferGF.hpp"
#include "SimdTransferGG.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ggk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ggk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 114263, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1215 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 114263, 85808, 6720, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1493, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1496, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1499, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1502, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1505, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1508, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1511, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1514, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1517, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1520, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1523, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1526, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1529, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1532, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1535, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1538, 3, 8, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1547, 3, 9, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1556, 3, 10, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1565, 3, 11, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1574, 3, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1583, 3, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1592, 3, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1601, 3, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1610, 3, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1619, 3, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1628, 3, 18, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1637, 3, 19, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1646, 3, 20, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1655, 3, 21, 62,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1664, 3, 23, 65,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1682, 3, 26, 71,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1700, 3, 29, 77,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1718, 3, 32, 83,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1736, 3, 35, 89,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1754, 3, 38, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1772, 3, 41, 101,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1790, 3, 44, 107,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1808, 3, 47, 113,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1826, 3, 50, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1844, 3, 53, 125,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1862, 3, 56, 131,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1880, 3, 59, 137,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1898, 3, 65, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1928, 3, 71, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1958, 3, 77, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1988, 3, 83, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2018, 3, 89, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2048, 3, 95, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2078, 3, 101, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2108, 3, 107, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2138, 3, 113, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2168, 3, 119, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2198, 3, 125, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2228, 3, 131, 253,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2258, 3, 143, 263,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2303, 3, 153, 278,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2348, 3, 163, 293,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2393, 3, 173, 308,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2438, 3, 183, 323,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2483, 3, 193, 338,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2528, 3, 203, 353,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2573, 3, 213, 368,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2618, 3, 223, 383,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2663, 3, 233, 398,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2708, 3, 243, 413,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2753, 3, 263, 428,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2816, 3, 278, 449,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2879, 3, 293, 470,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2942, 3, 308, 491,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3005, 3, 323, 512,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3068, 3, 338, 533,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3131, 3, 353, 554,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3194, 3, 368, 575,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3257, 3, 383, 596,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3320, 3, 398, 617,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3383, 3, 428, 638,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3467, 3, 449, 666,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3551, 3, 470, 694,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3635, 3, 491, 722,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3719, 3, 512, 750,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3803, 3, 533, 778,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3887, 3, 554, 806,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3971, 3, 575, 834,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4055, 3, 596, 862,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4139, 3, 638, 890,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4247, 3, 666, 926,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4355, 3, 694, 962,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4463, 3, 722, 998,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4571, 3, 750,
                                                                       1034, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4679, 3, 778,
                                                                       1070, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4787, 3, 806,
                                                                       1106, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4895, 3, 834,
                                                                       1142, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5003, 3, 890,
                                                                       1178, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5138, 3, 926,
                                                                       1223, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5273, 3, 962,
                                                                       1268, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5408, 3, 998,
                                                                       1313, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5543, 3, 1034,
                                                                       1358, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5678, 3, 1070,
                                                                       1403, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5813, 3, 1106,
                                                                       1448, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5948, 3, 8, 9,
                                                                       1499, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5954, 3, 9, 10,
                                                                       1502, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5960, 3, 10, 11,
                                                                       1505, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5966, 3, 11, 12,
                                                                       1508, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5972, 3, 12, 13,
                                                                       1511, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5978, 3, 13, 14,
                                                                       1514, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5984, 3, 14, 15,
                                                                       1517, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5990, 3, 15, 16,
                                                                       1520, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5996, 3, 16, 17,
                                                                       1523, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6002, 3, 17, 18,
                                                                       1526, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6008, 3, 18, 19,
                                                                       1529, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6014, 3, 19, 20,
                                                                       1532, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6020, 3, 20, 21,
                                                                       1535, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6026, 0, 3, 5948,
                                                                       1499, 5954, 1556, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6044, 0, 3, 5954,
                                                                       1502, 5960, 1565, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6062, 0, 3, 5960,
                                                                       1505, 5966, 1574, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6080, 0, 3, 5966,
                                                                       1508, 5972, 1583, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6098, 0, 3, 5972,
                                                                       1511, 5978, 1592, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6116, 0, 3, 5978,
                                                                       1514, 5984, 1601, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6134, 0, 3, 5984,
                                                                       1517, 5990, 1610, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6152, 0, 3, 5990,
                                                                       1520, 5996, 1619, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6170, 0, 3, 5996,
                                                                       1523, 6002, 1628, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6188, 0, 3, 6002,
                                                                       1526, 6008, 1637, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6206, 0, 3, 6008,
                                                                       1529, 6014, 1646, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6224, 0, 3, 6014,
                                                                       1532, 6020, 1655, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6242, 0, 3, 6026,
                                                                       1556, 6044, 65, 71, 1700,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6278, 0, 3, 6044,
                                                                       1565, 6062, 71, 77, 1718,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6314, 0, 3, 6062,
                                                                       1574, 6080, 77, 83, 1736,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6350, 0, 3, 6080,
                                                                       1583, 6098, 83, 89, 1754,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6386, 0, 3, 6098,
                                                                       1592, 6116, 89, 95, 1772,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6422, 0, 3, 6116,
                                                                       1601, 6134, 95, 101, 1790,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6458, 0, 3, 6134,
                                                                       1610, 6152, 101, 107,
                                                                       1808, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6494, 0, 3, 6152,
                                                                       1619, 6170, 107, 113,
                                                                       1826, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6530, 0, 3, 6170,
                                                                       1628, 6188, 113, 119,
                                                                       1844, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6566, 0, 3, 6188,
                                                                       1637, 6206, 119, 125,
                                                                       1862, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6602, 0, 3, 6206,
                                                                       1646, 6224, 125, 131,
                                                                       1880, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6638, 0, 3, 6242,
                                                                       1700, 6278, 143, 153,
                                                                       1958, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6698, 0, 3, 6278,
                                                                       1718, 6314, 153, 163,
                                                                       1988, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6758, 0, 3, 6314,
                                                                       1736, 6350, 163, 173,
                                                                       2018, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6818, 0, 3, 6350,
                                                                       1754, 6386, 173, 183,
                                                                       2048, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6878, 0, 3, 6386,
                                                                       1772, 6422, 183, 193,
                                                                       2078, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6938, 0, 3, 6422,
                                                                       1790, 6458, 193, 203,
                                                                       2108, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6998, 0, 3, 6458,
                                                                       1808, 6494, 203, 213,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7058, 0, 3, 6494,
                                                                       1826, 6530, 213, 223,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7118, 0, 3, 6530,
                                                                       1844, 6566, 223, 233,
                                                                       2198, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7178, 0, 3, 6566,
                                                                       1862, 6602, 233, 243,
                                                                       2228, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7238, 0, 3, 6638,
                                                                       1958, 6698, 263, 278,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7328, 0, 3, 6698,
                                                                       1988, 6758, 278, 293,
                                                                       2393, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7418, 0, 3, 6758,
                                                                       2018, 6818, 293, 308,
                                                                       2438, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7508, 0, 3, 6818,
                                                                       2048, 6878, 308, 323,
                                                                       2483, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7598, 0, 3, 6878,
                                                                       2078, 6938, 323, 338,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7688, 0, 3, 6938,
                                                                       2108, 6998, 338, 353,
                                                                       2573, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7778, 0, 3, 6998,
                                                                       2138, 7058, 353, 368,
                                                                       2618, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7868, 0, 3, 7058,
                                                                       2168, 7118, 368, 383,
                                                                       2663, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7958, 0, 3, 7118,
                                                                       2198, 7178, 383, 398,
                                                                       2708, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8048, 0, 3, 7238,
                                                                       2348, 7328, 428, 449,
                                                                       2879, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8174, 0, 3, 7328,
                                                                       2393, 7418, 449, 470,
                                                                       2942, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8300, 0, 3, 7418,
                                                                       2438, 7508, 470, 491,
                                                                       3005, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8426, 0, 3, 7508,
                                                                       2483, 7598, 491, 512,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8552, 0, 3, 7598,
                                                                       2528, 7688, 512, 533,
                                                                       3131, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8678, 0, 3, 7688,
                                                                       2573, 7778, 533, 554,
                                                                       3194, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8804, 0, 3, 7778,
                                                                       2618, 7868, 554, 575,
                                                                       3257, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8930, 0, 3, 7868,
                                                                       2663, 7958, 575, 596,
                                                                       3320, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9056, 0, 3, 8048,
                                                                       2879, 8174, 638, 666,
                                                                       3551, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9224, 0, 3, 8174,
                                                                       2942, 8300, 666, 694,
                                                                       3635, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9392, 0, 3, 8300,
                                                                       3005, 8426, 694, 722,
                                                                       3719, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9560, 0, 3, 8426,
                                                                       3068, 8552, 722, 750,
                                                                       3803, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9728, 0, 3, 8552,
                                                                       3131, 8678, 750, 778,
                                                                       3887, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9896, 0, 3, 8678,
                                                                       3194, 8804, 778, 806,
                                                                       3971, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10064, 0, 3, 8804,
                                                                       3257, 8930, 806, 834,
                                                                       4055, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10232, 0, 3, 9056,
                                                                       3551, 9224, 890, 926,
                                                                       4355, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10448, 0, 3, 9224,
                                                                       3635, 9392, 926, 962,
                                                                       4463, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10664, 0, 3, 9392,
                                                                       3719, 9560, 962, 998,
                                                                       4571, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10880, 0, 3, 9560,
                                                                       3803, 9728, 998, 1034,
                                                                       4679, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11096, 0, 3, 9728,
                                                                       3887, 9896, 1034, 1070,
                                                                       4787, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11312, 0, 3, 9896,
                                                                       3971, 10064, 1070, 1106,
                                                                       4895, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 11528, 0, 3,
                                                                       10232, 4355, 10448, 1178,
                                                                       1223, 5273, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 11798, 0, 3,
                                                                       10448, 4463, 10664, 1223,
                                                                       1268, 5408, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 12068, 0, 3,
                                                                       10664, 4571, 10880, 1268,
                                                                       1313, 5543, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 12338, 0, 3,
                                                                       10880, 4679, 11096, 1313,
                                                                       1358, 5678, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 12608, 0, 3,
                                                                       11096, 4787, 11312, 1358,
                                                                       1403, 5813, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12878, 3, 1493,
                                                                       1496, 5948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12888, 3, 1496,
                                                                       1499, 5954, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12898, 3, 1499,
                                                                       1502, 5960, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12908, 3, 1502,
                                                                       1505, 5966, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12918, 3, 1505,
                                                                       1508, 5972, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12928, 3, 1508,
                                                                       1511, 5978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12938, 3, 1511,
                                                                       1514, 5984, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12948, 3, 1514,
                                                                       1517, 5990, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12958, 3, 1517,
                                                                       1520, 5996, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12968, 3, 1520,
                                                                       1523, 6002, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12978, 3, 1523,
                                                                       1526, 6008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12988, 3, 1526,
                                                                       1529, 6014, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12998, 3, 1529,
                                                                       1532, 6020, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13008, 0, 3,
                                                                       12878, 5948, 12888, 1538,
                                                                       1547, 6026, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13038, 0, 3,
                                                                       12888, 5954, 12898, 1547,
                                                                       1556, 6044, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13068, 0, 3,
                                                                       12898, 5960, 12908, 1556,
                                                                       1565, 6062, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13098, 0, 3,
                                                                       12908, 5966, 12918, 1565,
                                                                       1574, 6080, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13128, 0, 3,
                                                                       12918, 5972, 12928, 1574,
                                                                       1583, 6098, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13158, 0, 3,
                                                                       12928, 5978, 12938, 1583,
                                                                       1592, 6116, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13188, 0, 3,
                                                                       12938, 5984, 12948, 1592,
                                                                       1601, 6134, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13218, 0, 3,
                                                                       12948, 5990, 12958, 1601,
                                                                       1610, 6152, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13248, 0, 3,
                                                                       12958, 5996, 12968, 1610,
                                                                       1619, 6170, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13278, 0, 3,
                                                                       12968, 6002, 12978, 1619,
                                                                       1628, 6188, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13308, 0, 3,
                                                                       12978, 6008, 12988, 1628,
                                                                       1637, 6206, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13338, 0, 3,
                                                                       12988, 6014, 12998, 1637,
                                                                       1646, 6224, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13368, 0, 3,
                                                                       13008, 6026, 13038, 1664,
                                                                       1682, 6242, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13428, 0, 3,
                                                                       13038, 6044, 13068, 1682,
                                                                       1700, 6278, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13488, 0, 3,
                                                                       13068, 6062, 13098, 1700,
                                                                       1718, 6314, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13548, 0, 3,
                                                                       13098, 6080, 13128, 1718,
                                                                       1736, 6350, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13608, 0, 3,
                                                                       13128, 6098, 13158, 1736,
                                                                       1754, 6386, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13668, 0, 3,
                                                                       13158, 6116, 13188, 1754,
                                                                       1772, 6422, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13728, 0, 3,
                                                                       13188, 6134, 13218, 1772,
                                                                       1790, 6458, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13788, 0, 3,
                                                                       13218, 6152, 13248, 1790,
                                                                       1808, 6494, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13848, 0, 3,
                                                                       13248, 6170, 13278, 1808,
                                                                       1826, 6530, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13908, 0, 3,
                                                                       13278, 6188, 13308, 1826,
                                                                       1844, 6566, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13968, 0, 3,
                                                                       13308, 6206, 13338, 1844,
                                                                       1862, 6602, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14028, 0, 3,
                                                                       13368, 6242, 13428, 1898,
                                                                       1928, 6638, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14128, 0, 3,
                                                                       13428, 6278, 13488, 1928,
                                                                       1958, 6698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14228, 0, 3,
                                                                       13488, 6314, 13548, 1958,
                                                                       1988, 6758, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14328, 0, 3,
                                                                       13548, 6350, 13608, 1988,
                                                                       2018, 6818, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14428, 0, 3,
                                                                       13608, 6386, 13668, 2018,
                                                                       2048, 6878, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14528, 0, 3,
                                                                       13668, 6422, 13728, 2048,
                                                                       2078, 6938, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14628, 0, 3,
                                                                       13728, 6458, 13788, 2078,
                                                                       2108, 6998, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14728, 0, 3,
                                                                       13788, 6494, 13848, 2108,
                                                                       2138, 7058, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14828, 0, 3,
                                                                       13848, 6530, 13908, 2138,
                                                                       2168, 7118, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14928, 0, 3,
                                                                       13908, 6566, 13968, 2168,
                                                                       2198, 7178, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15028, 0, 3,
                                                                       14028, 6638, 14128, 2258,
                                                                       2303, 7238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15178, 0, 3,
                                                                       14128, 6698, 14228, 2303,
                                                                       2348, 7328, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15328, 0, 3,
                                                                       14228, 6758, 14328, 2348,
                                                                       2393, 7418, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15478, 0, 3,
                                                                       14328, 6818, 14428, 2393,
                                                                       2438, 7508, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15628, 0, 3,
                                                                       14428, 6878, 14528, 2438,
                                                                       2483, 7598, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15778, 0, 3,
                                                                       14528, 6938, 14628, 2483,
                                                                       2528, 7688, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15928, 0, 3,
                                                                       14628, 6998, 14728, 2528,
                                                                       2573, 7778, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16078, 0, 3,
                                                                       14728, 7058, 14828, 2573,
                                                                       2618, 7868, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16228, 0, 3,
                                                                       14828, 7118, 14928, 2618,
                                                                       2663, 7958, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 16378, 0, 3,
                                                                       15028, 7238, 15178, 2753,
                                                                       2816, 8048, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 16588, 0, 3,
                                                                       15178, 7328, 15328, 2816,
                                                                       2879, 8174, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 16798, 0, 3,
                                                                       15328, 7418, 15478, 2879,
                                                                       2942, 8300, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 17008, 0, 3,
                                                                       15478, 7508, 15628, 2942,
                                                                       3005, 8426, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 17218, 0, 3,
                                                                       15628, 7598, 15778, 3005,
                                                                       3068, 8552, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 17428, 0, 3,
                                                                       15778, 7688, 15928, 3068,
                                                                       3131, 8678, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 17638, 0, 3,
                                                                       15928, 7778, 16078, 3131,
                                                                       3194, 8804, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 17848, 0, 3,
                                                                       16078, 7868, 16228, 3194,
                                                                       3257, 8930, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 18058, 0, 3,
                                                                       16378, 8048, 16588, 3383,
                                                                       3467, 9056, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 18338, 0, 3,
                                                                       16588, 8174, 16798, 3467,
                                                                       3551, 9224, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 18618, 0, 3,
                                                                       16798, 8300, 17008, 3551,
                                                                       3635, 9392, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 18898, 0, 3,
                                                                       17008, 8426, 17218, 3635,
                                                                       3719, 9560, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 19178, 0, 3,
                                                                       17218, 8552, 17428, 3719,
                                                                       3803, 9728, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 19458, 0, 3,
                                                                       17428, 8678, 17638, 3803,
                                                                       3887, 9896, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 19738, 0, 3,
                                                                       17638, 8804, 17848, 3887,
                                                                       3971, 10064, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 20018, 0, 3,
                                                                       18058, 9056, 18338, 4139,
                                                                       4247, 10232, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 20378, 0, 3,
                                                                       18338, 9224, 18618, 4247,
                                                                       4355, 10448, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 20738, 0, 3,
                                                                       18618, 9392, 18898, 4355,
                                                                       4463, 10664, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 21098, 0, 3,
                                                                       18898, 9560, 19178, 4463,
                                                                       4571, 10880, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 21458, 0, 3,
                                                                       19178, 9728, 19458, 4571,
                                                                       4679, 11096, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 21818, 0, 3,
                                                                       19458, 9896, 19738, 4679,
                                                                       4787, 11312, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 22178, 0, 3,
                                                                       20018, 10232, 20378, 5003,
                                                                       5138, 11528, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 22628, 0, 3,
                                                                       20378, 10448, 20738, 5138,
                                                                       5273, 11798, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 23078, 0, 3,
                                                                       20738, 10664, 21098, 5273,
                                                                       5408, 12068, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 23528, 0, 3,
                                                                       21098, 10880, 21458, 5408,
                                                                       5543, 12338, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 23978, 0, 3,
                                                                       21458, 11096, 21818, 5543,
                                                                       5678, 12608, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24428, 3, 5948,
                                                                       5954, 12898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24443, 3, 5954,
                                                                       5960, 12908, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24458, 3, 5960,
                                                                       5966, 12918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24473, 3, 5966,
                                                                       5972, 12928, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24488, 3, 5972,
                                                                       5978, 12938, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24503, 3, 5978,
                                                                       5984, 12948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24518, 3, 5984,
                                                                       5990, 12958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24533, 3, 5990,
                                                                       5996, 12968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24548, 3, 5996,
                                                                       6002, 12978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24563, 3, 6002,
                                                                       6008, 12988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24578, 3, 6008,
                                                                       6014, 12998, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24593, 0, 3,
                                                                       24428, 12898, 24443, 6026,
                                                                       6044, 13068, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24638, 0, 3,
                                                                       24443, 12908, 24458, 6044,
                                                                       6062, 13098, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24683, 0, 3,
                                                                       24458, 12918, 24473, 6062,
                                                                       6080, 13128, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24728, 0, 3,
                                                                       24473, 12928, 24488, 6080,
                                                                       6098, 13158, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24773, 0, 3,
                                                                       24488, 12938, 24503, 6098,
                                                                       6116, 13188, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24818, 0, 3,
                                                                       24503, 12948, 24518, 6116,
                                                                       6134, 13218, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24863, 0, 3,
                                                                       24518, 12958, 24533, 6134,
                                                                       6152, 13248, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24908, 0, 3,
                                                                       24533, 12968, 24548, 6152,
                                                                       6170, 13278, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24953, 0, 3,
                                                                       24548, 12978, 24563, 6170,
                                                                       6188, 13308, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24998, 0, 3,
                                                                       24563, 12988, 24578, 6188,
                                                                       6206, 13338, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25043, 0, 3,
                                                                       24593, 13068, 24638, 6242,
                                                                       6278, 13488, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25133, 0, 3,
                                                                       24638, 13098, 24683, 6278,
                                                                       6314, 13548, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25223, 0, 3,
                                                                       24683, 13128, 24728, 6314,
                                                                       6350, 13608, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25313, 0, 3,
                                                                       24728, 13158, 24773, 6350,
                                                                       6386, 13668, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25403, 0, 3,
                                                                       24773, 13188, 24818, 6386,
                                                                       6422, 13728, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25493, 0, 3,
                                                                       24818, 13218, 24863, 6422,
                                                                       6458, 13788, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25583, 0, 3,
                                                                       24863, 13248, 24908, 6458,
                                                                       6494, 13848, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25673, 0, 3,
                                                                       24908, 13278, 24953, 6494,
                                                                       6530, 13908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25763, 0, 3,
                                                                       24953, 13308, 24998, 6530,
                                                                       6566, 13968, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 25853, 0, 3,
                                                                       25043, 13488, 25133, 6638,
                                                                       6698, 14228, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26003, 0, 3,
                                                                       25133, 13548, 25223, 6698,
                                                                       6758, 14328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26153, 0, 3,
                                                                       25223, 13608, 25313, 6758,
                                                                       6818, 14428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26303, 0, 3,
                                                                       25313, 13668, 25403, 6818,
                                                                       6878, 14528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26453, 0, 3,
                                                                       25403, 13728, 25493, 6878,
                                                                       6938, 14628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26603, 0, 3,
                                                                       25493, 13788, 25583, 6938,
                                                                       6998, 14728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26753, 0, 3,
                                                                       25583, 13848, 25673, 6998,
                                                                       7058, 14828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26903, 0, 3,
                                                                       25673, 13908, 25763, 7058,
                                                                       7118, 14928, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27053, 0, 3,
                                                                       25853, 14228, 26003, 7238,
                                                                       7328, 15328, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27278, 0, 3,
                                                                       26003, 14328, 26153, 7328,
                                                                       7418, 15478, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27503, 0, 3,
                                                                       26153, 14428, 26303, 7418,
                                                                       7508, 15628, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27728, 0, 3,
                                                                       26303, 14528, 26453, 7508,
                                                                       7598, 15778, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27953, 0, 3,
                                                                       26453, 14628, 26603, 7598,
                                                                       7688, 15928, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28178, 0, 3,
                                                                       26603, 14728, 26753, 7688,
                                                                       7778, 16078, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28403, 0, 3,
                                                                       26753, 14828, 26903, 7778,
                                                                       7868, 16228, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 28628, 0, 3,
                                                                       27053, 15328, 27278, 8048,
                                                                       8174, 16798, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 28943, 0, 3,
                                                                       27278, 15478, 27503, 8174,
                                                                       8300, 17008, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 29258, 0, 3,
                                                                       27503, 15628, 27728, 8300,
                                                                       8426, 17218, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 29573, 0, 3,
                                                                       27728, 15778, 27953, 8426,
                                                                       8552, 17428, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 29888, 0, 3,
                                                                       27953, 15928, 28178, 8552,
                                                                       8678, 17638, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 30203, 0, 3,
                                                                       28178, 16078, 28403, 8678,
                                                                       8804, 17848, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 30518, 0, 3,
                                                                       28628, 16798, 28943, 9056,
                                                                       9224, 18618, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 30938, 0, 3,
                                                                       28943, 17008, 29258, 9224,
                                                                       9392, 18898, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 31358, 0, 3,
                                                                       29258, 17218, 29573, 9392,
                                                                       9560, 19178, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 31778, 0, 3,
                                                                       29573, 17428, 29888, 9560,
                                                                       9728, 19458, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 32198, 0, 3,
                                                                       29888, 17638, 30203, 9728,
                                                                       9896, 19738, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 32618, 0, 3,
                                                                       30518, 18618, 30938,
                                                                       10232, 10448, 20738,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 33158, 0, 3,
                                                                       30938, 18898, 31358,
                                                                       10448, 10664, 21098,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 33698, 0, 3,
                                                                       31358, 19178, 31778,
                                                                       10664, 10880, 21458,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 34238, 0, 3,
                                                                       31778, 19458, 32198,
                                                                       10880, 11096, 21818,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 34778, 0, 3,
                                                                       32618, 20738, 33158,
                                                                       11528, 11798, 23078,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 35453, 0, 3,
                                                                       33158, 21098, 33698,
                                                                       11798, 12068, 23528,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 36128, 0, 3,
                                                                       33698, 21458, 34238,
                                                                       12068, 12338, 23978,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36803, 3, 12878,
                                                                       12888, 24428, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36824, 3, 12888,
                                                                       12898, 24443, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36845, 3, 12898,
                                                                       12908, 24458, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36866, 3, 12908,
                                                                       12918, 24473, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36887, 3, 12918,
                                                                       12928, 24488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36908, 3, 12928,
                                                                       12938, 24503, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36929, 3, 12938,
                                                                       12948, 24518, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36950, 3, 12948,
                                                                       12958, 24533, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36971, 3, 12958,
                                                                       12968, 24548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36992, 3, 12968,
                                                                       12978, 24563, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 37013, 3, 12978,
                                                                       12988, 24578, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37034, 0, 3,
                                                                       36803, 24428, 36824,
                                                                       13008, 13038, 24593,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37097, 0, 3,
                                                                       36824, 24443, 36845,
                                                                       13038, 13068, 24638,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37160, 0, 3,
                                                                       36845, 24458, 36866,
                                                                       13068, 13098, 24683,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37223, 0, 3,
                                                                       36866, 24473, 36887,
                                                                       13098, 13128, 24728,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37286, 0, 3,
                                                                       36887, 24488, 36908,
                                                                       13128, 13158, 24773,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37349, 0, 3,
                                                                       36908, 24503, 36929,
                                                                       13158, 13188, 24818,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37412, 0, 3,
                                                                       36929, 24518, 36950,
                                                                       13188, 13218, 24863,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37475, 0, 3,
                                                                       36950, 24533, 36971,
                                                                       13218, 13248, 24908,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37538, 0, 3,
                                                                       36971, 24548, 36992,
                                                                       13248, 13278, 24953,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37601, 0, 3,
                                                                       36992, 24563, 37013,
                                                                       13278, 13308, 24998,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37664, 0, 3,
                                                                       37034, 24593, 37097,
                                                                       13368, 13428, 25043,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37790, 0, 3,
                                                                       37097, 24638, 37160,
                                                                       13428, 13488, 25133,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37916, 0, 3,
                                                                       37160, 24683, 37223,
                                                                       13488, 13548, 25223,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38042, 0, 3,
                                                                       37223, 24728, 37286,
                                                                       13548, 13608, 25313,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38168, 0, 3,
                                                                       37286, 24773, 37349,
                                                                       13608, 13668, 25403,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38294, 0, 3,
                                                                       37349, 24818, 37412,
                                                                       13668, 13728, 25493,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38420, 0, 3,
                                                                       37412, 24863, 37475,
                                                                       13728, 13788, 25583,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38546, 0, 3,
                                                                       37475, 24908, 37538,
                                                                       13788, 13848, 25673,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38672, 0, 3,
                                                                       37538, 24953, 37601,
                                                                       13848, 13908, 25763,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 38798, 0, 3,
                                                                       37664, 25043, 37790,
                                                                       14028, 14128, 25853,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39008, 0, 3,
                                                                       37790, 25133, 37916,
                                                                       14128, 14228, 26003,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39218, 0, 3,
                                                                       37916, 25223, 38042,
                                                                       14228, 14328, 26153,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39428, 0, 3,
                                                                       38042, 25313, 38168,
                                                                       14328, 14428, 26303,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39638, 0, 3,
                                                                       38168, 25403, 38294,
                                                                       14428, 14528, 26453,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39848, 0, 3,
                                                                       38294, 25493, 38420,
                                                                       14528, 14628, 26603,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 40058, 0, 3,
                                                                       38420, 25583, 38546,
                                                                       14628, 14728, 26753,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 40268, 0, 3,
                                                                       38546, 25673, 38672,
                                                                       14728, 14828, 26903,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 40478, 0, 3,
                                                                       38798, 25853, 39008,
                                                                       15028, 15178, 27053,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 40793, 0, 3,
                                                                       39008, 26003, 39218,
                                                                       15178, 15328, 27278,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41108, 0, 3,
                                                                       39218, 26153, 39428,
                                                                       15328, 15478, 27503,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41423, 0, 3,
                                                                       39428, 26303, 39638,
                                                                       15478, 15628, 27728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41738, 0, 3,
                                                                       39638, 26453, 39848,
                                                                       15628, 15778, 27953,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 42053, 0, 3,
                                                                       39848, 26603, 40058,
                                                                       15778, 15928, 28178,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 42368, 0, 3,
                                                                       40058, 26753, 40268,
                                                                       15928, 16078, 28403,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 42683, 0, 3,
                                                                       40478, 27053, 40793,
                                                                       16378, 16588, 28628,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 43124, 0, 3,
                                                                       40793, 27278, 41108,
                                                                       16588, 16798, 28943,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 43565, 0, 3,
                                                                       41108, 27503, 41423,
                                                                       16798, 17008, 29258,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 44006, 0, 3,
                                                                       41423, 27728, 41738,
                                                                       17008, 17218, 29573,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 44447, 0, 3,
                                                                       41738, 27953, 42053,
                                                                       17218, 17428, 29888,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 44888, 0, 3,
                                                                       42053, 28178, 42368,
                                                                       17428, 17638, 30203,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 45329, 0, 3,
                                                                       42683, 28628, 43124,
                                                                       18058, 18338, 30518,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 45917, 0, 3,
                                                                       43124, 28943, 43565,
                                                                       18338, 18618, 30938,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 46505, 0, 3,
                                                                       43565, 29258, 44006,
                                                                       18618, 18898, 31358,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 47093, 0, 3,
                                                                       44006, 29573, 44447,
                                                                       18898, 19178, 31778,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 47681, 0, 3,
                                                                       44447, 29888, 44888,
                                                                       19178, 19458, 32198,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 48269, 0, 3,
                                                                       45329, 30518, 45917,
                                                                       20018, 20378, 32618,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 49025, 0, 3,
                                                                       45917, 30938, 46505,
                                                                       20378, 20738, 33158,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 49781, 0, 3,
                                                                       46505, 31358, 47093,
                                                                       20738, 21098, 33698,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 50537, 0, 3,
                                                                       47093, 31778, 47681,
                                                                       21098, 21458, 34238,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 51293, 0, 3,
                                                                       48269, 32618, 49025,
                                                                       22178, 22628, 34778,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 52238, 0, 3,
                                                                       49025, 33158, 49781,
                                                                       22628, 23078, 35453,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 53183, 0, 3,
                                                                       49781, 33698, 50537,
                                                                       23078, 23528, 36128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54128, 3, 24428,
                                                                       24443, 36845, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54156, 3, 24443,
                                                                       24458, 36866, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54184, 3, 24458,
                                                                       24473, 36887, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54212, 3, 24473,
                                                                       24488, 36908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54240, 3, 24488,
                                                                       24503, 36929, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54268, 3, 24503,
                                                                       24518, 36950, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54296, 3, 24518,
                                                                       24533, 36971, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54324, 3, 24533,
                                                                       24548, 36992, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54352, 3, 24548,
                                                                       24563, 37013, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54380, 0, 3,
                                                                       54128, 36845, 54156,
                                                                       24593, 24638, 37160,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54464, 0, 3,
                                                                       54156, 36866, 54184,
                                                                       24638, 24683, 37223,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54548, 0, 3,
                                                                       54184, 36887, 54212,
                                                                       24683, 24728, 37286,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54632, 0, 3,
                                                                       54212, 36908, 54240,
                                                                       24728, 24773, 37349,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54716, 0, 3,
                                                                       54240, 36929, 54268,
                                                                       24773, 24818, 37412,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54800, 0, 3,
                                                                       54268, 36950, 54296,
                                                                       24818, 24863, 37475,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54884, 0, 3,
                                                                       54296, 36971, 54324,
                                                                       24863, 24908, 37538,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54968, 0, 3,
                                                                       54324, 36992, 54352,
                                                                       24908, 24953, 37601,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55052, 0, 3,
                                                                       54380, 37160, 54464,
                                                                       25043, 25133, 37916,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55220, 0, 3,
                                                                       54464, 37223, 54548,
                                                                       25133, 25223, 38042,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55388, 0, 3,
                                                                       54548, 37286, 54632,
                                                                       25223, 25313, 38168,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55556, 0, 3,
                                                                       54632, 37349, 54716,
                                                                       25313, 25403, 38294,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55724, 0, 3,
                                                                       54716, 37412, 54800,
                                                                       25403, 25493, 38420,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55892, 0, 3,
                                                                       54800, 37475, 54884,
                                                                       25493, 25583, 38546,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 56060, 0, 3,
                                                                       54884, 37538, 54968,
                                                                       25583, 25673, 38672,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 56228, 0, 3,
                                                                       55052, 37916, 55220,
                                                                       25853, 26003, 39218,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 56508, 0, 3,
                                                                       55220, 38042, 55388,
                                                                       26003, 26153, 39428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 56788, 0, 3,
                                                                       55388, 38168, 55556,
                                                                       26153, 26303, 39638,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 57068, 0, 3,
                                                                       55556, 38294, 55724,
                                                                       26303, 26453, 39848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 57348, 0, 3,
                                                                       55724, 38420, 55892,
                                                                       26453, 26603, 40058,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 57628, 0, 3,
                                                                       55892, 38546, 56060,
                                                                       26603, 26753, 40268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 57908, 0, 3,
                                                                       56228, 39218, 56508,
                                                                       27053, 27278, 41108,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 58328, 0, 3,
                                                                       56508, 39428, 56788,
                                                                       27278, 27503, 41423,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 58748, 0, 3,
                                                                       56788, 39638, 57068,
                                                                       27503, 27728, 41738,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 59168, 0, 3,
                                                                       57068, 39848, 57348,
                                                                       27728, 27953, 42053,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 59588, 0, 3,
                                                                       57348, 40058, 57628,
                                                                       27953, 28178, 42368,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 60008, 0, 3,
                                                                       57908, 41108, 58328,
                                                                       28628, 28943, 43565,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 60596, 0, 3,
                                                                       58328, 41423, 58748,
                                                                       28943, 29258, 44006,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 61184, 0, 3,
                                                                       58748, 41738, 59168,
                                                                       29258, 29573, 44447,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 61772, 0, 3,
                                                                       59168, 42053, 59588,
                                                                       29573, 29888, 44888,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 62360, 0, 3,
                                                                       60008, 43565, 60596,
                                                                       30518, 30938, 46505,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 63144, 0, 3,
                                                                       60596, 44006, 61184,
                                                                       30938, 31358, 47093,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 63928, 0, 3,
                                                                       61184, 44447, 61772,
                                                                       31358, 31778, 47681,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 64712, 0, 3,
                                                                       62360, 46505, 63144,
                                                                       32618, 33158, 49781,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 65720, 0, 3,
                                                                       63144, 47093, 63928,
                                                                       33158, 33698, 50537,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 66728, 0, 3,
                                                                       64712, 49781, 65720,
                                                                       34778, 35453, 53183,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 67988, 3, 36803,
                                                                       36824, 54128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68024, 3, 36824,
                                                                       36845, 54156, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68060, 3, 36845,
                                                                       36866, 54184, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68096, 3, 36866,
                                                                       36887, 54212, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68132, 3, 36887,
                                                                       36908, 54240, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68168, 3, 36908,
                                                                       36929, 54268, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68204, 3, 36929,
                                                                       36950, 54296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68240, 3, 36950,
                                                                       36971, 54324, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68276, 3, 36971,
                                                                       36992, 54352, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68312, 0, 3,
                                                                       67988, 54128, 68024,
                                                                       37034, 37097, 54380,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68420, 0, 3,
                                                                       68024, 54156, 68060,
                                                                       37097, 37160, 54464,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68528, 0, 3,
                                                                       68060, 54184, 68096,
                                                                       37160, 37223, 54548,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68636, 0, 3,
                                                                       68096, 54212, 68132,
                                                                       37223, 37286, 54632,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68744, 0, 3,
                                                                       68132, 54240, 68168,
                                                                       37286, 37349, 54716,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68852, 0, 3,
                                                                       68168, 54268, 68204,
                                                                       37349, 37412, 54800,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68960, 0, 3,
                                                                       68204, 54296, 68240,
                                                                       37412, 37475, 54884,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 69068, 0, 3,
                                                                       68240, 54324, 68276,
                                                                       37475, 37538, 54968,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 69176, 0, 3,
                                                                       68312, 54380, 68420,
                                                                       37664, 37790, 55052,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 69392, 0, 3,
                                                                       68420, 54464, 68528,
                                                                       37790, 37916, 55220,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 69608, 0, 3,
                                                                       68528, 54548, 68636,
                                                                       37916, 38042, 55388,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 69824, 0, 3,
                                                                       68636, 54632, 68744,
                                                                       38042, 38168, 55556,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 70040, 0, 3,
                                                                       68744, 54716, 68852,
                                                                       38168, 38294, 55724,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 70256, 0, 3,
                                                                       68852, 54800, 68960,
                                                                       38294, 38420, 55892,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 70472, 0, 3,
                                                                       68960, 54884, 69068,
                                                                       38420, 38546, 56060,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 70688, 0, 3,
                                                                       69176, 55052, 69392,
                                                                       38798, 39008, 56228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 71048, 0, 3,
                                                                       69392, 55220, 69608,
                                                                       39008, 39218, 56508,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 71408, 0, 3,
                                                                       69608, 55388, 69824,
                                                                       39218, 39428, 56788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 71768, 0, 3,
                                                                       69824, 55556, 70040,
                                                                       39428, 39638, 57068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 72128, 0, 3,
                                                                       70040, 55724, 70256,
                                                                       39638, 39848, 57348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 72488, 0, 3,
                                                                       70256, 55892, 70472,
                                                                       39848, 40058, 57628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 72848, 0, 3,
                                                                       70688, 56228, 71048,
                                                                       40478, 40793, 57908,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 73388, 0, 3,
                                                                       71048, 56508, 71408,
                                                                       40793, 41108, 58328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 73928, 0, 3,
                                                                       71408, 56788, 71768,
                                                                       41108, 41423, 58748,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 74468, 0, 3,
                                                                       71768, 57068, 72128,
                                                                       41423, 41738, 59168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 75008, 0, 3,
                                                                       72128, 57348, 72488,
                                                                       41738, 42053, 59588,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 75548, 0, 3,
                                                                       72848, 57908, 73388,
                                                                       42683, 43124, 60008,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 76304, 0, 3,
                                                                       73388, 58328, 73928,
                                                                       43124, 43565, 60596,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 77060, 0, 3,
                                                                       73928, 58748, 74468,
                                                                       43565, 44006, 61184,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 77816, 0, 3,
                                                                       74468, 59168, 75008,
                                                                       44006, 44447, 61772,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 78572, 0, 3,
                                                                       75548, 60008, 76304,
                                                                       45329, 45917, 62360,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 79580, 0, 3,
                                                                       76304, 60596, 77060,
                                                                       45917, 46505, 63144,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 80588, 0, 3,
                                                                       77060, 61184, 77816,
                                                                       46505, 47093, 63928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 81596, 0, 3,
                                                                       78572, 62360, 79580,
                                                                       48269, 49025, 64712,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 82892, 0, 3,
                                                                       79580, 63144, 80588,
                                                                       49025, 49781, 65720,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 84188, 0, 3,
                                                                       81596, 64712, 82892,
                                                                       51293, 52238, 66728,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 85808, 72848, 540, ncols);

                    simdfunc::contract_primitives(buffer, 86573, 75548, 756, ncols);

                    simdfunc::contract_primitives(buffer, 87644, 78572, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 89072, 81596, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 90908, 84188, 1620, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 86348, 85808, 15, 1, nmax);

        simdtrf::transform_k_inner(buffer, 87329, 86573, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 88652, 87644, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 90368, 89072, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 92528, 90908, 45, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 93203, 86348, 87329, 15,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 93878, 87329, 88652, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 94823, 88652, 90368, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 96083, 90368, 92528, 15,
                                             nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 97703, 93203, 93878, 15,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 99053, 93878, 94823, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 100943, 94823, 96083, 15,
                                             nmax);

        simdtrf::compute_hrr_gf_out_of_first(buffer, coordinates, 103463, 97703, 99053, 15,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 105713, 99053, 100943, 15,
                                             nmax);

        simdtrf::compute_hrr_gg_out_of_first(buffer, coordinates, 108863, 103463, 105713, 15,
                                             nmax);

        simdtrf::transform_g_inner(buffer, 112238, 108863, 15, 15, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 112238, 135, nmax);
    }

    for (size_t m = 0; m < 1215; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
