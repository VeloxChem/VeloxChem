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


#include "SimdThreeCenterElectronRepulsionRecHFI.hpp"

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
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hfi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hfi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 73876, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 73876, 56768, 4745, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 14,
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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1493, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1496, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1499, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1502, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1505, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1508, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1511, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1514, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1517, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1520, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1523, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1526, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1529, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1532, 3, 10, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1541, 3, 11, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1550, 3, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1559, 3, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1568, 3, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1577, 3, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1586, 3, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1595, 3, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1604, 3, 18, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1613, 3, 19, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1622, 3, 20, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1631, 3, 21, 62,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1640, 3, 29, 77,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1658, 3, 32, 83,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1676, 3, 35, 89,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1694, 3, 38, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1712, 3, 41, 101,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1730, 3, 44, 107,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1748, 3, 47, 113,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1766, 3, 50, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1784, 3, 53, 125,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1802, 3, 56, 131,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1820, 3, 59, 137,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1838, 3, 77, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1868, 3, 83, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1898, 3, 89, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1928, 3, 95, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1958, 3, 101, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1988, 3, 107, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2018, 3, 113, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2048, 3, 119, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2078, 3, 125, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2108, 3, 131, 253,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2138, 3, 163, 293,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2183, 3, 173, 308,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2228, 3, 183, 323,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2273, 3, 193, 338,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2318, 3, 203, 353,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2363, 3, 213, 368,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2408, 3, 223, 383,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2453, 3, 233, 398,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2498, 3, 243, 413,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2543, 3, 293, 470,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2606, 3, 308, 491,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2669, 3, 323, 512,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2732, 3, 338, 533,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2795, 3, 353, 554,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2858, 3, 368, 575,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2921, 3, 383, 596,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2984, 3, 398, 617,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3047, 3, 470, 694,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3131, 3, 491, 722,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3215, 3, 512, 750,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3299, 3, 533, 778,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3383, 3, 554, 806,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3467, 3, 575, 834,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3551, 3, 596, 862,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3635, 3, 694, 962,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3743, 3, 722, 998,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3851, 3, 750,
                                                                       1034, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3959, 3, 778,
                                                                       1070, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4067, 3, 806,
                                                                       1106, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4175, 3, 834,
                                                                       1142, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4283, 3, 962,
                                                                       1268, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4418, 3, 998,
                                                                       1313, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4553, 3, 1034,
                                                                       1358, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4688, 3, 1070,
                                                                       1403, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4823, 3, 1106,
                                                                       1448, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4958, 3, 8, 9,
                                                                       1493, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4964, 3, 9, 10,
                                                                       1496, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4970, 3, 10, 11,
                                                                       1499, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4976, 3, 11, 12,
                                                                       1502, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4982, 3, 12, 13,
                                                                       1505, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4988, 3, 13, 14,
                                                                       1508, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4994, 3, 14, 15,
                                                                       1511, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5000, 3, 15, 16,
                                                                       1514, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5006, 3, 16, 17,
                                                                       1517, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5012, 3, 17, 18,
                                                                       1520, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5018, 3, 18, 19,
                                                                       1523, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5024, 3, 19, 20,
                                                                       1526, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5030, 3, 20, 21,
                                                                       1529, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5036, 0, 3, 4958,
                                                                       1493, 4964, 1532, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5054, 0, 3, 4964,
                                                                       1496, 4970, 1541, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5072, 0, 3, 4970,
                                                                       1499, 4976, 1550, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5090, 0, 3, 4976,
                                                                       1502, 4982, 1559, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5108, 0, 3, 4982,
                                                                       1505, 4988, 1568, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5126, 0, 3, 4988,
                                                                       1508, 4994, 1577, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5144, 0, 3, 4994,
                                                                       1511, 5000, 1586, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5162, 0, 3, 5000,
                                                                       1514, 5006, 1595, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5180, 0, 3, 5006,
                                                                       1517, 5012, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5198, 0, 3, 5012,
                                                                       1520, 5018, 1613, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5216, 0, 3, 5018,
                                                                       1523, 5024, 1622, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5234, 0, 3, 5024,
                                                                       1526, 5030, 1631, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5252, 0, 3, 5036,
                                                                       1532, 5054, 65, 71, 1640,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 5054,
                                                                       1541, 5072, 71, 77, 1658,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5324, 0, 3, 5072,
                                                                       1550, 5090, 77, 83, 1676,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5360, 0, 3, 5090,
                                                                       1559, 5108, 83, 89, 1694,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5396, 0, 3, 5108,
                                                                       1568, 5126, 89, 95, 1712,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5432, 0, 3, 5126,
                                                                       1577, 5144, 95, 101, 1730,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5468, 0, 3, 5144,
                                                                       1586, 5162, 101, 107,
                                                                       1748, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5504, 0, 3, 5162,
                                                                       1595, 5180, 107, 113,
                                                                       1766, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5540, 0, 3, 5180,
                                                                       1604, 5198, 113, 119,
                                                                       1784, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5576, 0, 3, 5198,
                                                                       1613, 5216, 119, 125,
                                                                       1802, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5612, 0, 3, 5216,
                                                                       1622, 5234, 125, 131,
                                                                       1820, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5648, 0, 3, 5252,
                                                                       1640, 5288, 143, 153,
                                                                       1838, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5708, 0, 3, 5288,
                                                                       1658, 5324, 153, 163,
                                                                       1868, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5768, 0, 3, 5324,
                                                                       1676, 5360, 163, 173,
                                                                       1898, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5828, 0, 3, 5360,
                                                                       1694, 5396, 173, 183,
                                                                       1928, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5888, 0, 3, 5396,
                                                                       1712, 5432, 183, 193,
                                                                       1958, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5948, 0, 3, 5432,
                                                                       1730, 5468, 193, 203,
                                                                       1988, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6008, 0, 3, 5468,
                                                                       1748, 5504, 203, 213,
                                                                       2018, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6068, 0, 3, 5504,
                                                                       1766, 5540, 213, 223,
                                                                       2048, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6128, 0, 3, 5540,
                                                                       1784, 5576, 223, 233,
                                                                       2078, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6188, 0, 3, 5576,
                                                                       1802, 5612, 233, 243,
                                                                       2108, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 5648,
                                                                       1838, 5708, 263, 278,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6338, 0, 3, 5708,
                                                                       1868, 5768, 278, 293,
                                                                       2183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6428, 0, 3, 5768,
                                                                       1898, 5828, 293, 308,
                                                                       2228, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6518, 0, 3, 5828,
                                                                       1928, 5888, 308, 323,
                                                                       2273, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 5888,
                                                                       1958, 5948, 323, 338,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6698, 0, 3, 5948,
                                                                       1988, 6008, 338, 353,
                                                                       2363, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6788, 0, 3, 6008,
                                                                       2018, 6068, 353, 368,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6878, 0, 3, 6068,
                                                                       2048, 6128, 368, 383,
                                                                       2453, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6968, 0, 3, 6128,
                                                                       2078, 6188, 383, 398,
                                                                       2498, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7058, 0, 3, 6248,
                                                                       2138, 6338, 428, 449,
                                                                       2543, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7184, 0, 3, 6338,
                                                                       2183, 6428, 449, 470,
                                                                       2606, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7310, 0, 3, 6428,
                                                                       2228, 6518, 470, 491,
                                                                       2669, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7436, 0, 3, 6518,
                                                                       2273, 6608, 491, 512,
                                                                       2732, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7562, 0, 3, 6608,
                                                                       2318, 6698, 512, 533,
                                                                       2795, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7688, 0, 3, 6698,
                                                                       2363, 6788, 533, 554,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7814, 0, 3, 6788,
                                                                       2408, 6878, 554, 575,
                                                                       2921, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7940, 0, 3, 6878,
                                                                       2453, 6968, 575, 596,
                                                                       2984, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8066, 0, 3, 7058,
                                                                       2543, 7184, 638, 666,
                                                                       3047, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8234, 0, 3, 7184,
                                                                       2606, 7310, 666, 694,
                                                                       3131, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8402, 0, 3, 7310,
                                                                       2669, 7436, 694, 722,
                                                                       3215, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8570, 0, 3, 7436,
                                                                       2732, 7562, 722, 750,
                                                                       3299, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8738, 0, 3, 7562,
                                                                       2795, 7688, 750, 778,
                                                                       3383, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8906, 0, 3, 7688,
                                                                       2858, 7814, 778, 806,
                                                                       3467, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9074, 0, 3, 7814,
                                                                       2921, 7940, 806, 834,
                                                                       3551, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9242, 0, 3, 8066,
                                                                       3047, 8234, 890, 926,
                                                                       3635, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9458, 0, 3, 8234,
                                                                       3131, 8402, 926, 962,
                                                                       3743, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9674, 0, 3, 8402,
                                                                       3215, 8570, 962, 998,
                                                                       3851, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9890, 0, 3, 8570,
                                                                       3299, 8738, 998, 1034,
                                                                       3959, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10106, 0, 3, 8738,
                                                                       3383, 8906, 1034, 1070,
                                                                       4067, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10322, 0, 3, 8906,
                                                                       3467, 9074, 1070, 1106,
                                                                       4175, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 10538, 0, 3, 9242,
                                                                       3635, 9458, 1178, 1223,
                                                                       4283, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 10808, 0, 3, 9458,
                                                                       3743, 9674, 1223, 1268,
                                                                       4418, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 11078, 0, 3, 9674,
                                                                       3851, 9890, 1268, 1313,
                                                                       4553, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 11348, 0, 3, 9890,
                                                                       3959, 10106, 1313, 1358,
                                                                       4688, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 11618, 0, 3,
                                                                       10106, 4067, 10322, 1358,
                                                                       1403, 4823, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11888, 3, 1493,
                                                                       1496, 4970, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11898, 3, 1496,
                                                                       1499, 4976, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11908, 3, 1499,
                                                                       1502, 4982, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11918, 3, 1502,
                                                                       1505, 4988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11928, 3, 1505,
                                                                       1508, 4994, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11938, 3, 1508,
                                                                       1511, 5000, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11948, 3, 1511,
                                                                       1514, 5006, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11958, 3, 1514,
                                                                       1517, 5012, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11968, 3, 1517,
                                                                       1520, 5018, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11978, 3, 1520,
                                                                       1523, 5024, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11988, 3, 1523,
                                                                       1526, 5030, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 11998, 0, 3,
                                                                       11888, 4970, 11898, 1532,
                                                                       1541, 5072, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12028, 0, 3,
                                                                       11898, 4976, 11908, 1541,
                                                                       1550, 5090, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12058, 0, 3,
                                                                       11908, 4982, 11918, 1550,
                                                                       1559, 5108, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12088, 0, 3,
                                                                       11918, 4988, 11928, 1559,
                                                                       1568, 5126, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12118, 0, 3,
                                                                       11928, 4994, 11938, 1568,
                                                                       1577, 5144, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12148, 0, 3,
                                                                       11938, 5000, 11948, 1577,
                                                                       1586, 5162, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12178, 0, 3,
                                                                       11948, 5006, 11958, 1586,
                                                                       1595, 5180, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12208, 0, 3,
                                                                       11958, 5012, 11968, 1595,
                                                                       1604, 5198, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12238, 0, 3,
                                                                       11968, 5018, 11978, 1604,
                                                                       1613, 5216, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12268, 0, 3,
                                                                       11978, 5024, 11988, 1613,
                                                                       1622, 5234, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12298, 0, 3,
                                                                       11998, 5072, 12028, 1640,
                                                                       1658, 5324, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12358, 0, 3,
                                                                       12028, 5090, 12058, 1658,
                                                                       1676, 5360, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12418, 0, 3,
                                                                       12058, 5108, 12088, 1676,
                                                                       1694, 5396, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12478, 0, 3,
                                                                       12088, 5126, 12118, 1694,
                                                                       1712, 5432, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12538, 0, 3,
                                                                       12118, 5144, 12148, 1712,
                                                                       1730, 5468, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12598, 0, 3,
                                                                       12148, 5162, 12178, 1730,
                                                                       1748, 5504, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12658, 0, 3,
                                                                       12178, 5180, 12208, 1748,
                                                                       1766, 5540, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12718, 0, 3,
                                                                       12208, 5198, 12238, 1766,
                                                                       1784, 5576, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12778, 0, 3,
                                                                       12238, 5216, 12268, 1784,
                                                                       1802, 5612, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 12838, 0, 3,
                                                                       12298, 5324, 12358, 1838,
                                                                       1868, 5768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 12938, 0, 3,
                                                                       12358, 5360, 12418, 1868,
                                                                       1898, 5828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13038, 0, 3,
                                                                       12418, 5396, 12478, 1898,
                                                                       1928, 5888, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13138, 0, 3,
                                                                       12478, 5432, 12538, 1928,
                                                                       1958, 5948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13238, 0, 3,
                                                                       12538, 5468, 12598, 1958,
                                                                       1988, 6008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13338, 0, 3,
                                                                       12598, 5504, 12658, 1988,
                                                                       2018, 6068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13438, 0, 3,
                                                                       12658, 5540, 12718, 2018,
                                                                       2048, 6128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13538, 0, 3,
                                                                       12718, 5576, 12778, 2048,
                                                                       2078, 6188, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13638, 0, 3,
                                                                       12838, 5768, 12938, 2138,
                                                                       2183, 6428, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13788, 0, 3,
                                                                       12938, 5828, 13038, 2183,
                                                                       2228, 6518, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13938, 0, 3,
                                                                       13038, 5888, 13138, 2228,
                                                                       2273, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 14088, 0, 3,
                                                                       13138, 5948, 13238, 2273,
                                                                       2318, 6698, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 14238, 0, 3,
                                                                       13238, 6008, 13338, 2318,
                                                                       2363, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 14388, 0, 3,
                                                                       13338, 6068, 13438, 2363,
                                                                       2408, 6878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 14538, 0, 3,
                                                                       13438, 6128, 13538, 2408,
                                                                       2453, 6968, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 14688, 0, 3,
                                                                       13638, 6428, 13788, 2543,
                                                                       2606, 7310, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 14898, 0, 3,
                                                                       13788, 6518, 13938, 2606,
                                                                       2669, 7436, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 15108, 0, 3,
                                                                       13938, 6608, 14088, 2669,
                                                                       2732, 7562, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 15318, 0, 3,
                                                                       14088, 6698, 14238, 2732,
                                                                       2795, 7688, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 15528, 0, 3,
                                                                       14238, 6788, 14388, 2795,
                                                                       2858, 7814, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 15738, 0, 3,
                                                                       14388, 6878, 14538, 2858,
                                                                       2921, 7940, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 15948, 0, 3,
                                                                       14688, 7310, 14898, 3047,
                                                                       3131, 8402, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 16228, 0, 3,
                                                                       14898, 7436, 15108, 3131,
                                                                       3215, 8570, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 16508, 0, 3,
                                                                       15108, 7562, 15318, 3215,
                                                                       3299, 8738, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 16788, 0, 3,
                                                                       15318, 7688, 15528, 3299,
                                                                       3383, 8906, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 17068, 0, 3,
                                                                       15528, 7814, 15738, 3383,
                                                                       3467, 9074, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 17348, 0, 3,
                                                                       15948, 8402, 16228, 3635,
                                                                       3743, 9674, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 17708, 0, 3,
                                                                       16228, 8570, 16508, 3743,
                                                                       3851, 9890, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 18068, 0, 3,
                                                                       16508, 8738, 16788, 3851,
                                                                       3959, 10106, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 18428, 0, 3,
                                                                       16788, 8906, 17068, 3959,
                                                                       4067, 10322, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 18788, 0, 3,
                                                                       17348, 9674, 17708, 4283,
                                                                       4418, 11078, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 19238, 0, 3,
                                                                       17708, 9890, 18068, 4418,
                                                                       4553, 11348, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 19688, 0, 3,
                                                                       18068, 10106, 18428, 4553,
                                                                       4688, 11618, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20138, 3, 4958,
                                                                       4964, 11888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20153, 3, 4964,
                                                                       4970, 11898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20168, 3, 4970,
                                                                       4976, 11908, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20183, 3, 4976,
                                                                       4982, 11918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20198, 3, 4982,
                                                                       4988, 11928, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20213, 3, 4988,
                                                                       4994, 11938, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20228, 3, 4994,
                                                                       5000, 11948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20243, 3, 5000,
                                                                       5006, 11958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20258, 3, 5006,
                                                                       5012, 11968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20273, 3, 5012,
                                                                       5018, 11978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20288, 3, 5018,
                                                                       5024, 11988, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20303, 0, 3,
                                                                       20138, 11888, 20153, 5036,
                                                                       5054, 11998, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20348, 0, 3,
                                                                       20153, 11898, 20168, 5054,
                                                                       5072, 12028, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20393, 0, 3,
                                                                       20168, 11908, 20183, 5072,
                                                                       5090, 12058, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20438, 0, 3,
                                                                       20183, 11918, 20198, 5090,
                                                                       5108, 12088, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20483, 0, 3,
                                                                       20198, 11928, 20213, 5108,
                                                                       5126, 12118, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20528, 0, 3,
                                                                       20213, 11938, 20228, 5126,
                                                                       5144, 12148, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20573, 0, 3,
                                                                       20228, 11948, 20243, 5144,
                                                                       5162, 12178, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20618, 0, 3,
                                                                       20243, 11958, 20258, 5162,
                                                                       5180, 12208, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20663, 0, 3,
                                                                       20258, 11968, 20273, 5180,
                                                                       5198, 12238, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20708, 0, 3,
                                                                       20273, 11978, 20288, 5198,
                                                                       5216, 12268, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 20753, 0, 3,
                                                                       20303, 11998, 20348, 5252,
                                                                       5288, 12298, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 20843, 0, 3,
                                                                       20348, 12028, 20393, 5288,
                                                                       5324, 12358, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 20933, 0, 3,
                                                                       20393, 12058, 20438, 5324,
                                                                       5360, 12418, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21023, 0, 3,
                                                                       20438, 12088, 20483, 5360,
                                                                       5396, 12478, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21113, 0, 3,
                                                                       20483, 12118, 20528, 5396,
                                                                       5432, 12538, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21203, 0, 3,
                                                                       20528, 12148, 20573, 5432,
                                                                       5468, 12598, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21293, 0, 3,
                                                                       20573, 12178, 20618, 5468,
                                                                       5504, 12658, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21383, 0, 3,
                                                                       20618, 12208, 20663, 5504,
                                                                       5540, 12718, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21473, 0, 3,
                                                                       20663, 12238, 20708, 5540,
                                                                       5576, 12778, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 21563, 0, 3,
                                                                       20753, 12298, 20843, 5648,
                                                                       5708, 12838, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 21713, 0, 3,
                                                                       20843, 12358, 20933, 5708,
                                                                       5768, 12938, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 21863, 0, 3,
                                                                       20933, 12418, 21023, 5768,
                                                                       5828, 13038, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22013, 0, 3,
                                                                       21023, 12478, 21113, 5828,
                                                                       5888, 13138, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22163, 0, 3,
                                                                       21113, 12538, 21203, 5888,
                                                                       5948, 13238, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22313, 0, 3,
                                                                       21203, 12598, 21293, 5948,
                                                                       6008, 13338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22463, 0, 3,
                                                                       21293, 12658, 21383, 6008,
                                                                       6068, 13438, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22613, 0, 3,
                                                                       21383, 12718, 21473, 6068,
                                                                       6128, 13538, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 22763, 0, 3,
                                                                       21563, 12838, 21713, 6248,
                                                                       6338, 13638, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 22988, 0, 3,
                                                                       21713, 12938, 21863, 6338,
                                                                       6428, 13788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 23213, 0, 3,
                                                                       21863, 13038, 22013, 6428,
                                                                       6518, 13938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 23438, 0, 3,
                                                                       22013, 13138, 22163, 6518,
                                                                       6608, 14088, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 23663, 0, 3,
                                                                       22163, 13238, 22313, 6608,
                                                                       6698, 14238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 23888, 0, 3,
                                                                       22313, 13338, 22463, 6698,
                                                                       6788, 14388, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 24113, 0, 3,
                                                                       22463, 13438, 22613, 6788,
                                                                       6878, 14538, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 24338, 0, 3,
                                                                       22763, 13638, 22988, 7058,
                                                                       7184, 14688, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 24653, 0, 3,
                                                                       22988, 13788, 23213, 7184,
                                                                       7310, 14898, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 24968, 0, 3,
                                                                       23213, 13938, 23438, 7310,
                                                                       7436, 15108, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 25283, 0, 3,
                                                                       23438, 14088, 23663, 7436,
                                                                       7562, 15318, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 25598, 0, 3,
                                                                       23663, 14238, 23888, 7562,
                                                                       7688, 15528, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 25913, 0, 3,
                                                                       23888, 14388, 24113, 7688,
                                                                       7814, 15738, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 26228, 0, 3,
                                                                       24338, 14688, 24653, 8066,
                                                                       8234, 15948, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 26648, 0, 3,
                                                                       24653, 14898, 24968, 8234,
                                                                       8402, 16228, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 27068, 0, 3,
                                                                       24968, 15108, 25283, 8402,
                                                                       8570, 16508, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 27488, 0, 3,
                                                                       25283, 15318, 25598, 8570,
                                                                       8738, 16788, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 27908, 0, 3,
                                                                       25598, 15528, 25913, 8738,
                                                                       8906, 17068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 28328, 0, 3,
                                                                       26228, 15948, 26648, 9242,
                                                                       9458, 17348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 28868, 0, 3,
                                                                       26648, 16228, 27068, 9458,
                                                                       9674, 17708, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 29408, 0, 3,
                                                                       27068, 16508, 27488, 9674,
                                                                       9890, 18068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 29948, 0, 3,
                                                                       27488, 16788, 27908, 9890,
                                                                       10106, 18428, ncols,
                                                                       gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 30488, 0, 3,
                                                                       28328, 17348, 28868,
                                                                       10538, 10808, 18788,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 31163, 0, 3,
                                                                       28868, 17708, 29408,
                                                                       10808, 11078, 19238,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 31838, 0, 3,
                                                                       29408, 18068, 29948,
                                                                       11078, 11348, 19688,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32513, 3, 11888,
                                                                       11898, 20168, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32534, 3, 11898,
                                                                       11908, 20183, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32555, 3, 11908,
                                                                       11918, 20198, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32576, 3, 11918,
                                                                       11928, 20213, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32597, 3, 11928,
                                                                       11938, 20228, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32618, 3, 11938,
                                                                       11948, 20243, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32639, 3, 11948,
                                                                       11958, 20258, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32660, 3, 11958,
                                                                       11968, 20273, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32681, 3, 11968,
                                                                       11978, 20288, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 32702, 0, 3,
                                                                       32513, 20168, 32534,
                                                                       11998, 12028, 20393,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 32765, 0, 3,
                                                                       32534, 20183, 32555,
                                                                       12028, 12058, 20438,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 32828, 0, 3,
                                                                       32555, 20198, 32576,
                                                                       12058, 12088, 20483,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 32891, 0, 3,
                                                                       32576, 20213, 32597,
                                                                       12088, 12118, 20528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 32954, 0, 3,
                                                                       32597, 20228, 32618,
                                                                       12118, 12148, 20573,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 33017, 0, 3,
                                                                       32618, 20243, 32639,
                                                                       12148, 12178, 20618,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 33080, 0, 3,
                                                                       32639, 20258, 32660,
                                                                       12178, 12208, 20663,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 33143, 0, 3,
                                                                       32660, 20273, 32681,
                                                                       12208, 12238, 20708,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33206, 0, 3,
                                                                       32702, 20393, 32765,
                                                                       12298, 12358, 20933,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33332, 0, 3,
                                                                       32765, 20438, 32828,
                                                                       12358, 12418, 21023,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33458, 0, 3,
                                                                       32828, 20483, 32891,
                                                                       12418, 12478, 21113,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33584, 0, 3,
                                                                       32891, 20528, 32954,
                                                                       12478, 12538, 21203,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33710, 0, 3,
                                                                       32954, 20573, 33017,
                                                                       12538, 12598, 21293,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33836, 0, 3,
                                                                       33017, 20618, 33080,
                                                                       12598, 12658, 21383,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33962, 0, 3,
                                                                       33080, 20663, 33143,
                                                                       12658, 12718, 21473,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 34088, 0, 3,
                                                                       33206, 20933, 33332,
                                                                       12838, 12938, 21863,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 34298, 0, 3,
                                                                       33332, 21023, 33458,
                                                                       12938, 13038, 22013,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 34508, 0, 3,
                                                                       33458, 21113, 33584,
                                                                       13038, 13138, 22163,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 34718, 0, 3,
                                                                       33584, 21203, 33710,
                                                                       13138, 13238, 22313,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 34928, 0, 3,
                                                                       33710, 21293, 33836,
                                                                       13238, 13338, 22463,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 35138, 0, 3,
                                                                       33836, 21383, 33962,
                                                                       13338, 13438, 22613,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 35348, 0, 3,
                                                                       34088, 21863, 34298,
                                                                       13638, 13788, 23213,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 35663, 0, 3,
                                                                       34298, 22013, 34508,
                                                                       13788, 13938, 23438,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 35978, 0, 3,
                                                                       34508, 22163, 34718,
                                                                       13938, 14088, 23663,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 36293, 0, 3,
                                                                       34718, 22313, 34928,
                                                                       14088, 14238, 23888,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 36608, 0, 3,
                                                                       34928, 22463, 35138,
                                                                       14238, 14388, 24113,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 36923, 0, 3,
                                                                       35348, 23213, 35663,
                                                                       14688, 14898, 24968,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 37364, 0, 3,
                                                                       35663, 23438, 35978,
                                                                       14898, 15108, 25283,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 37805, 0, 3,
                                                                       35978, 23663, 36293,
                                                                       15108, 15318, 25598,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 38246, 0, 3,
                                                                       36293, 23888, 36608,
                                                                       15318, 15528, 25913,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 38687, 0, 3,
                                                                       36923, 24968, 37364,
                                                                       15948, 16228, 27068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 39275, 0, 3,
                                                                       37364, 25283, 37805,
                                                                       16228, 16508, 27488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 39863, 0, 3,
                                                                       37805, 25598, 38246,
                                                                       16508, 16788, 27908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 40451, 0, 3,
                                                                       38687, 27068, 39275,
                                                                       17348, 17708, 29408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 41207, 0, 3,
                                                                       39275, 27488, 39863,
                                                                       17708, 18068, 29948,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 41963, 0, 3,
                                                                       40451, 29408, 41207,
                                                                       18788, 19238, 31838,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42908, 3, 20138,
                                                                       20153, 32513, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42936, 3, 20153,
                                                                       20168, 32534, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42964, 3, 20168,
                                                                       20183, 32555, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42992, 3, 20183,
                                                                       20198, 32576, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43020, 3, 20198,
                                                                       20213, 32597, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43048, 3, 20213,
                                                                       20228, 32618, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43076, 3, 20228,
                                                                       20243, 32639, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43104, 3, 20243,
                                                                       20258, 32660, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43132, 3, 20258,
                                                                       20273, 32681, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43160, 0, 3,
                                                                       42908, 32513, 42936,
                                                                       20303, 20348, 32702,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43244, 0, 3,
                                                                       42936, 32534, 42964,
                                                                       20348, 20393, 32765,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43328, 0, 3,
                                                                       42964, 32555, 42992,
                                                                       20393, 20438, 32828,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43412, 0, 3,
                                                                       42992, 32576, 43020,
                                                                       20438, 20483, 32891,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43496, 0, 3,
                                                                       43020, 32597, 43048,
                                                                       20483, 20528, 32954,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43580, 0, 3,
                                                                       43048, 32618, 43076,
                                                                       20528, 20573, 33017,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43664, 0, 3,
                                                                       43076, 32639, 43104,
                                                                       20573, 20618, 33080,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43748, 0, 3,
                                                                       43104, 32660, 43132,
                                                                       20618, 20663, 33143,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 43832, 0, 3,
                                                                       43160, 32702, 43244,
                                                                       20753, 20843, 33206,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 44000, 0, 3,
                                                                       43244, 32765, 43328,
                                                                       20843, 20933, 33332,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 44168, 0, 3,
                                                                       43328, 32828, 43412,
                                                                       20933, 21023, 33458,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 44336, 0, 3,
                                                                       43412, 32891, 43496,
                                                                       21023, 21113, 33584,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 44504, 0, 3,
                                                                       43496, 32954, 43580,
                                                                       21113, 21203, 33710,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 44672, 0, 3,
                                                                       43580, 33017, 43664,
                                                                       21203, 21293, 33836,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 44840, 0, 3,
                                                                       43664, 33080, 43748,
                                                                       21293, 21383, 33962,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 45008, 0, 3,
                                                                       43832, 33206, 44000,
                                                                       21563, 21713, 34088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 45288, 0, 3,
                                                                       44000, 33332, 44168,
                                                                       21713, 21863, 34298,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 45568, 0, 3,
                                                                       44168, 33458, 44336,
                                                                       21863, 22013, 34508,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 45848, 0, 3,
                                                                       44336, 33584, 44504,
                                                                       22013, 22163, 34718,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 46128, 0, 3,
                                                                       44504, 33710, 44672,
                                                                       22163, 22313, 34928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 46408, 0, 3,
                                                                       44672, 33836, 44840,
                                                                       22313, 22463, 35138,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 46688, 0, 3,
                                                                       45008, 34088, 45288,
                                                                       22763, 22988, 35348,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 47108, 0, 3,
                                                                       45288, 34298, 45568,
                                                                       22988, 23213, 35663,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 47528, 0, 3,
                                                                       45568, 34508, 45848,
                                                                       23213, 23438, 35978,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 47948, 0, 3,
                                                                       45848, 34718, 46128,
                                                                       23438, 23663, 36293,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 48368, 0, 3,
                                                                       46128, 34928, 46408,
                                                                       23663, 23888, 36608,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 48788, 0, 3,
                                                                       46688, 35348, 47108,
                                                                       24338, 24653, 36923,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 49376, 0, 3,
                                                                       47108, 35663, 47528,
                                                                       24653, 24968, 37364,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 49964, 0, 3,
                                                                       47528, 35978, 47948,
                                                                       24968, 25283, 37805,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 50552, 0, 3,
                                                                       47948, 36293, 48368,
                                                                       25283, 25598, 38246,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 51140, 0, 3,
                                                                       48788, 36923, 49376,
                                                                       26228, 26648, 38687,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 51924, 0, 3,
                                                                       49376, 37364, 49964,
                                                                       26648, 27068, 39275,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 52708, 0, 3,
                                                                       49964, 37805, 50552,
                                                                       27068, 27488, 39863,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 53492, 0, 3,
                                                                       51140, 38687, 51924,
                                                                       28328, 28868, 40451,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 54500, 0, 3,
                                                                       51924, 39275, 52708,
                                                                       28868, 29408, 41207,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 55508, 0, 3,
                                                                       53492, 40451, 54500,
                                                                       30488, 31163, 41963,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 56768, 48788, 588, ncols);

                    simdfunc::contract_primitives(buffer, 57629, 51140, 784, ncols);

                    simdfunc::contract_primitives(buffer, 58777, 53492, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 60253, 55508, 1260, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 57356, 56768, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 58413, 57629, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 59785, 58777, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 61513, 60253, 45, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 62098, 57356, 58413, 13,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 62917, 58413, 59785, 13,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 64009, 59785, 61513, 13,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 65413, 62098, 62917, 13,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 67051, 62917, 64009, 13,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 69235, 65413, 67051, 13,
                                             nmax);

        simdtrf::transform_f_inner(buffer, 71965, 69235, 21, 13, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 71965, 91, nmax);
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
