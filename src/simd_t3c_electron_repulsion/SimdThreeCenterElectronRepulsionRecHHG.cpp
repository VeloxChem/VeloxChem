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


#include "SimdThreeCenterElectronRepulsionRecHHG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hhg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hhg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 79687, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1089 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 79687, 41335, 5430, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2153, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2156, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2159, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2162, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2165, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2168, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2171, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2174, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2177, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2180, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2183, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2186, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2189, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2192, 3, 10, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2201, 3, 11, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2210, 3, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2219, 3, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2228, 3, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2237, 3, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2246, 3, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2255, 3, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2264, 3, 18, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2273, 3, 19, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2282, 3, 20, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2291, 3, 21, 62,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2300, 3, 29, 77,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2318, 3, 32, 83,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2336, 3, 35, 89,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2354, 3, 38, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2372, 3, 41, 101,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2390, 3, 44, 107,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2408, 3, 47, 113,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2426, 3, 50, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2444, 3, 53, 125,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2462, 3, 56, 131,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2480, 3, 59, 137,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2498, 3, 77, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2528, 3, 83, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2558, 3, 89, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2588, 3, 95, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2618, 3, 101, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2648, 3, 107, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2678, 3, 113, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2708, 3, 119, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2738, 3, 125, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2768, 3, 131, 253,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2798, 3, 163, 293,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2843, 3, 173, 308,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2888, 3, 183, 323,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2933, 3, 193, 338,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2978, 3, 203, 353,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3023, 3, 213, 368,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3068, 3, 223, 383,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3113, 3, 233, 398,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3158, 3, 243, 413,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3203, 3, 293, 470,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3266, 3, 308, 491,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3329, 3, 323, 512,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3392, 3, 338, 533,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3455, 3, 353, 554,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3518, 3, 368, 575,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3581, 3, 383, 596,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3644, 3, 398, 617,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3707, 3, 470, 694,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3791, 3, 491, 722,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3875, 3, 512, 750,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3959, 3, 533, 778,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4043, 3, 554, 806,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4127, 3, 575, 834,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4211, 3, 596, 862,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4295, 3, 694, 962,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4403, 3, 722, 998,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4511, 3, 750,
                                                                       1034, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4619, 3, 778,
                                                                       1070, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4727, 3, 806,
                                                                       1106, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4835, 3, 834,
                                                                       1142, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4943, 3, 962,
                                                                       1268, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5078, 3, 998,
                                                                       1313, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5213, 3, 1034,
                                                                       1358, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5348, 3, 1070,
                                                                       1403, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5483, 3, 1106,
                                                                       1448, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 5618, 3, 1268,
                                                                       1603, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 5783, 3, 1313,
                                                                       1658, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 5948, 3, 1358,
                                                                       1713, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6113, 3, 1403,
                                                                       1768, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 6278, 3, 1603,
                                                                       1955, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 6476, 3, 1658,
                                                                       2021, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 6674, 3, 1713,
                                                                       2087, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6872, 3, 8, 9,
                                                                       2153, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6878, 3, 9, 10,
                                                                       2156, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6884, 3, 10, 11,
                                                                       2159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6890, 3, 11, 12,
                                                                       2162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6896, 3, 12, 13,
                                                                       2165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6902, 3, 13, 14,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6908, 3, 14, 15,
                                                                       2171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6914, 3, 15, 16,
                                                                       2174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6920, 3, 16, 17,
                                                                       2177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6926, 3, 17, 18,
                                                                       2180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6932, 3, 18, 19,
                                                                       2183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6938, 3, 19, 20,
                                                                       2186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6944, 3, 20, 21,
                                                                       2189, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6950, 0, 3, 6872,
                                                                       2153, 6878, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6968, 0, 3, 6878,
                                                                       2156, 6884, 2201, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6986, 0, 3, 6884,
                                                                       2159, 6890, 2210, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7004, 0, 3, 6890,
                                                                       2162, 6896, 2219, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7022, 0, 3, 6896,
                                                                       2165, 6902, 2228, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7040, 0, 3, 6902,
                                                                       2168, 6908, 2237, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7058, 0, 3, 6908,
                                                                       2171, 6914, 2246, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7076, 0, 3, 6914,
                                                                       2174, 6920, 2255, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7094, 0, 3, 6920,
                                                                       2177, 6926, 2264, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7112, 0, 3, 6926,
                                                                       2180, 6932, 2273, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7130, 0, 3, 6932,
                                                                       2183, 6938, 2282, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 6938,
                                                                       2186, 6944, 2291, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7166, 0, 3, 6950,
                                                                       2192, 6968, 65, 71, 2300,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7202, 0, 3, 6968,
                                                                       2201, 6986, 71, 77, 2318,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7238, 0, 3, 6986,
                                                                       2210, 7004, 77, 83, 2336,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7274, 0, 3, 7004,
                                                                       2219, 7022, 83, 89, 2354,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7310, 0, 3, 7022,
                                                                       2228, 7040, 89, 95, 2372,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7346, 0, 3, 7040,
                                                                       2237, 7058, 95, 101, 2390,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7382, 0, 3, 7058,
                                                                       2246, 7076, 101, 107,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7418, 0, 3, 7076,
                                                                       2255, 7094, 107, 113,
                                                                       2426, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7454, 0, 3, 7094,
                                                                       2264, 7112, 113, 119,
                                                                       2444, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7490, 0, 3, 7112,
                                                                       2273, 7130, 119, 125,
                                                                       2462, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7526, 0, 3, 7130,
                                                                       2282, 7148, 125, 131,
                                                                       2480, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7562, 0, 3, 7166,
                                                                       2300, 7202, 143, 153,
                                                                       2498, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7622, 0, 3, 7202,
                                                                       2318, 7238, 153, 163,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7682, 0, 3, 7238,
                                                                       2336, 7274, 163, 173,
                                                                       2558, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7742, 0, 3, 7274,
                                                                       2354, 7310, 173, 183,
                                                                       2588, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7802, 0, 3, 7310,
                                                                       2372, 7346, 183, 193,
                                                                       2618, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7862, 0, 3, 7346,
                                                                       2390, 7382, 193, 203,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7922, 0, 3, 7382,
                                                                       2408, 7418, 203, 213,
                                                                       2678, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7982, 0, 3, 7418,
                                                                       2426, 7454, 213, 223,
                                                                       2708, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8042, 0, 3, 7454,
                                                                       2444, 7490, 223, 233,
                                                                       2738, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8102, 0, 3, 7490,
                                                                       2462, 7526, 233, 243,
                                                                       2768, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8162, 0, 3, 7562,
                                                                       2498, 7622, 263, 278,
                                                                       2798, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8252, 0, 3, 7622,
                                                                       2528, 7682, 278, 293,
                                                                       2843, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8342, 0, 3, 7682,
                                                                       2558, 7742, 293, 308,
                                                                       2888, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8432, 0, 3, 7742,
                                                                       2588, 7802, 308, 323,
                                                                       2933, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8522, 0, 3, 7802,
                                                                       2618, 7862, 323, 338,
                                                                       2978, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8612, 0, 3, 7862,
                                                                       2648, 7922, 338, 353,
                                                                       3023, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8702, 0, 3, 7922,
                                                                       2678, 7982, 353, 368,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8792, 0, 3, 7982,
                                                                       2708, 8042, 368, 383,
                                                                       3113, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8882, 0, 3, 8042,
                                                                       2738, 8102, 383, 398,
                                                                       3158, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8972, 0, 3, 8162,
                                                                       2798, 8252, 428, 449,
                                                                       3203, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9098, 0, 3, 8252,
                                                                       2843, 8342, 449, 470,
                                                                       3266, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9224, 0, 3, 8342,
                                                                       2888, 8432, 470, 491,
                                                                       3329, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9350, 0, 3, 8432,
                                                                       2933, 8522, 491, 512,
                                                                       3392, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9476, 0, 3, 8522,
                                                                       2978, 8612, 512, 533,
                                                                       3455, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9602, 0, 3, 8612,
                                                                       3023, 8702, 533, 554,
                                                                       3518, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9728, 0, 3, 8702,
                                                                       3068, 8792, 554, 575,
                                                                       3581, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9854, 0, 3, 8792,
                                                                       3113, 8882, 575, 596,
                                                                       3644, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9980, 0, 3, 8972,
                                                                       3203, 9098, 638, 666,
                                                                       3707, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10148, 0, 3, 9098,
                                                                       3266, 9224, 666, 694,
                                                                       3791, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10316, 0, 3, 9224,
                                                                       3329, 9350, 694, 722,
                                                                       3875, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10484, 0, 3, 9350,
                                                                       3392, 9476, 722, 750,
                                                                       3959, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10652, 0, 3, 9476,
                                                                       3455, 9602, 750, 778,
                                                                       4043, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10820, 0, 3, 9602,
                                                                       3518, 9728, 778, 806,
                                                                       4127, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10988, 0, 3, 9728,
                                                                       3581, 9854, 806, 834,
                                                                       4211, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11156, 0, 3, 9980,
                                                                       3707, 10148, 890, 926,
                                                                       4295, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11372, 0, 3,
                                                                       10148, 3791, 10316, 926,
                                                                       962, 4403, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11588, 0, 3,
                                                                       10316, 3875, 10484, 962,
                                                                       998, 4511, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11804, 0, 3,
                                                                       10484, 3959, 10652, 998,
                                                                       1034, 4619, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12020, 0, 3,
                                                                       10652, 4043, 10820, 1034,
                                                                       1070, 4727, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12236, 0, 3,
                                                                       10820, 4127, 10988, 1070,
                                                                       1106, 4835, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 12452, 0, 3,
                                                                       11156, 4295, 11372, 1178,
                                                                       1223, 4943, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 12722, 0, 3,
                                                                       11372, 4403, 11588, 1223,
                                                                       1268, 5078, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 12992, 0, 3,
                                                                       11588, 4511, 11804, 1268,
                                                                       1313, 5213, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13262, 0, 3,
                                                                       11804, 4619, 12020, 1313,
                                                                       1358, 5348, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13532, 0, 3,
                                                                       12020, 4727, 12236, 1358,
                                                                       1403, 5483, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 13802, 0, 3,
                                                                       12452, 4943, 12722, 1493,
                                                                       1548, 5618, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 14132, 0, 3,
                                                                       12722, 5078, 12992, 1548,
                                                                       1603, 5783, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 14462, 0, 3,
                                                                       12992, 5213, 13262, 1603,
                                                                       1658, 5948, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 14792, 0, 3,
                                                                       13262, 5348, 13532, 1658,
                                                                       1713, 6113, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 15122, 0, 3,
                                                                       13802, 5618, 14132, 1823,
                                                                       1889, 6278, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 15518, 0, 3,
                                                                       14132, 5783, 14462, 1889,
                                                                       1955, 6476, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 15914, 0, 3,
                                                                       14462, 5948, 14792, 1955,
                                                                       2021, 6674, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16310, 3, 2153,
                                                                       2156, 6884, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16320, 3, 2156,
                                                                       2159, 6890, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16330, 3, 2159,
                                                                       2162, 6896, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16340, 3, 2162,
                                                                       2165, 6902, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16350, 3, 2165,
                                                                       2168, 6908, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16360, 3, 2168,
                                                                       2171, 6914, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16370, 3, 2171,
                                                                       2174, 6920, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16380, 3, 2174,
                                                                       2177, 6926, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16390, 3, 2177,
                                                                       2180, 6932, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16400, 3, 2180,
                                                                       2183, 6938, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16410, 3, 2183,
                                                                       2186, 6944, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16420, 0, 3,
                                                                       16310, 6884, 16320, 2192,
                                                                       2201, 6986, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16450, 0, 3,
                                                                       16320, 6890, 16330, 2201,
                                                                       2210, 7004, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16480, 0, 3,
                                                                       16330, 6896, 16340, 2210,
                                                                       2219, 7022, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16510, 0, 3,
                                                                       16340, 6902, 16350, 2219,
                                                                       2228, 7040, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16540, 0, 3,
                                                                       16350, 6908, 16360, 2228,
                                                                       2237, 7058, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16570, 0, 3,
                                                                       16360, 6914, 16370, 2237,
                                                                       2246, 7076, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16600, 0, 3,
                                                                       16370, 6920, 16380, 2246,
                                                                       2255, 7094, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16630, 0, 3,
                                                                       16380, 6926, 16390, 2255,
                                                                       2264, 7112, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16660, 0, 3,
                                                                       16390, 6932, 16400, 2264,
                                                                       2273, 7130, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16690, 0, 3,
                                                                       16400, 6938, 16410, 2273,
                                                                       2282, 7148, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16720, 0, 3,
                                                                       16420, 6986, 16450, 2300,
                                                                       2318, 7238, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16780, 0, 3,
                                                                       16450, 7004, 16480, 2318,
                                                                       2336, 7274, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16840, 0, 3,
                                                                       16480, 7022, 16510, 2336,
                                                                       2354, 7310, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16900, 0, 3,
                                                                       16510, 7040, 16540, 2354,
                                                                       2372, 7346, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16960, 0, 3,
                                                                       16540, 7058, 16570, 2372,
                                                                       2390, 7382, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17020, 0, 3,
                                                                       16570, 7076, 16600, 2390,
                                                                       2408, 7418, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17080, 0, 3,
                                                                       16600, 7094, 16630, 2408,
                                                                       2426, 7454, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17140, 0, 3,
                                                                       16630, 7112, 16660, 2426,
                                                                       2444, 7490, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17200, 0, 3,
                                                                       16660, 7130, 16690, 2444,
                                                                       2462, 7526, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17260, 0, 3,
                                                                       16720, 7238, 16780, 2498,
                                                                       2528, 7682, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17360, 0, 3,
                                                                       16780, 7274, 16840, 2528,
                                                                       2558, 7742, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17460, 0, 3,
                                                                       16840, 7310, 16900, 2558,
                                                                       2588, 7802, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17560, 0, 3,
                                                                       16900, 7346, 16960, 2588,
                                                                       2618, 7862, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17660, 0, 3,
                                                                       16960, 7382, 17020, 2618,
                                                                       2648, 7922, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17760, 0, 3,
                                                                       17020, 7418, 17080, 2648,
                                                                       2678, 7982, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17860, 0, 3,
                                                                       17080, 7454, 17140, 2678,
                                                                       2708, 8042, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17960, 0, 3,
                                                                       17140, 7490, 17200, 2708,
                                                                       2738, 8102, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18060, 0, 3,
                                                                       17260, 7682, 17360, 2798,
                                                                       2843, 8342, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18210, 0, 3,
                                                                       17360, 7742, 17460, 2843,
                                                                       2888, 8432, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18360, 0, 3,
                                                                       17460, 7802, 17560, 2888,
                                                                       2933, 8522, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18510, 0, 3,
                                                                       17560, 7862, 17660, 2933,
                                                                       2978, 8612, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18660, 0, 3,
                                                                       17660, 7922, 17760, 2978,
                                                                       3023, 8702, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18810, 0, 3,
                                                                       17760, 7982, 17860, 3023,
                                                                       3068, 8792, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18960, 0, 3,
                                                                       17860, 8042, 17960, 3068,
                                                                       3113, 8882, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19110, 0, 3,
                                                                       18060, 8342, 18210, 3203,
                                                                       3266, 9224, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19320, 0, 3,
                                                                       18210, 8432, 18360, 3266,
                                                                       3329, 9350, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19530, 0, 3,
                                                                       18360, 8522, 18510, 3329,
                                                                       3392, 9476, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19740, 0, 3,
                                                                       18510, 8612, 18660, 3392,
                                                                       3455, 9602, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19950, 0, 3,
                                                                       18660, 8702, 18810, 3455,
                                                                       3518, 9728, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20160, 0, 3,
                                                                       18810, 8792, 18960, 3518,
                                                                       3581, 9854, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20370, 0, 3,
                                                                       19110, 9224, 19320, 3707,
                                                                       3791, 10316, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20650, 0, 3,
                                                                       19320, 9350, 19530, 3791,
                                                                       3875, 10484, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20930, 0, 3,
                                                                       19530, 9476, 19740, 3875,
                                                                       3959, 10652, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21210, 0, 3,
                                                                       19740, 9602, 19950, 3959,
                                                                       4043, 10820, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21490, 0, 3,
                                                                       19950, 9728, 20160, 4043,
                                                                       4127, 10988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 21770, 0, 3,
                                                                       20370, 10316, 20650, 4295,
                                                                       4403, 11588, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 22130, 0, 3,
                                                                       20650, 10484, 20930, 4403,
                                                                       4511, 11804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 22490, 0, 3,
                                                                       20930, 10652, 21210, 4511,
                                                                       4619, 12020, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 22850, 0, 3,
                                                                       21210, 10820, 21490, 4619,
                                                                       4727, 12236, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 23210, 0, 3,
                                                                       21770, 11588, 22130, 4943,
                                                                       5078, 12992, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 23660, 0, 3,
                                                                       22130, 11804, 22490, 5078,
                                                                       5213, 13262, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 24110, 0, 3,
                                                                       22490, 12020, 22850, 5213,
                                                                       5348, 13532, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 24560, 0, 3,
                                                                       23210, 12992, 23660, 5618,
                                                                       5783, 14462, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 25110, 0, 3,
                                                                       23660, 13262, 24110, 5783,
                                                                       5948, 14792, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 25660, 0, 3,
                                                                       24560, 14462, 25110, 6278,
                                                                       6476, 15914, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26320, 3, 6872,
                                                                       6878, 16310, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26335, 3, 6878,
                                                                       6884, 16320, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26350, 3, 6884,
                                                                       6890, 16330, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26365, 3, 6890,
                                                                       6896, 16340, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26380, 3, 6896,
                                                                       6902, 16350, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26395, 3, 6902,
                                                                       6908, 16360, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26410, 3, 6908,
                                                                       6914, 16370, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26425, 3, 6914,
                                                                       6920, 16380, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26440, 3, 6920,
                                                                       6926, 16390, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26455, 3, 6926,
                                                                       6932, 16400, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26470, 3, 6932,
                                                                       6938, 16410, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26485, 0, 3,
                                                                       26320, 16310, 26335, 6950,
                                                                       6968, 16420, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26530, 0, 3,
                                                                       26335, 16320, 26350, 6968,
                                                                       6986, 16450, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26575, 0, 3,
                                                                       26350, 16330, 26365, 6986,
                                                                       7004, 16480, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26620, 0, 3,
                                                                       26365, 16340, 26380, 7004,
                                                                       7022, 16510, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26665, 0, 3,
                                                                       26380, 16350, 26395, 7022,
                                                                       7040, 16540, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26710, 0, 3,
                                                                       26395, 16360, 26410, 7040,
                                                                       7058, 16570, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26755, 0, 3,
                                                                       26410, 16370, 26425, 7058,
                                                                       7076, 16600, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26800, 0, 3,
                                                                       26425, 16380, 26440, 7076,
                                                                       7094, 16630, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26845, 0, 3,
                                                                       26440, 16390, 26455, 7094,
                                                                       7112, 16660, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26890, 0, 3,
                                                                       26455, 16400, 26470, 7112,
                                                                       7130, 16690, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 26935, 0, 3,
                                                                       26485, 16420, 26530, 7166,
                                                                       7202, 16720, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27025, 0, 3,
                                                                       26530, 16450, 26575, 7202,
                                                                       7238, 16780, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27115, 0, 3,
                                                                       26575, 16480, 26620, 7238,
                                                                       7274, 16840, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27205, 0, 3,
                                                                       26620, 16510, 26665, 7274,
                                                                       7310, 16900, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27295, 0, 3,
                                                                       26665, 16540, 26710, 7310,
                                                                       7346, 16960, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27385, 0, 3,
                                                                       26710, 16570, 26755, 7346,
                                                                       7382, 17020, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27475, 0, 3,
                                                                       26755, 16600, 26800, 7382,
                                                                       7418, 17080, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27565, 0, 3,
                                                                       26800, 16630, 26845, 7418,
                                                                       7454, 17140, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27655, 0, 3,
                                                                       26845, 16660, 26890, 7454,
                                                                       7490, 17200, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 27745, 0, 3,
                                                                       26935, 16720, 27025, 7562,
                                                                       7622, 17260, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 27895, 0, 3,
                                                                       27025, 16780, 27115, 7622,
                                                                       7682, 17360, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28045, 0, 3,
                                                                       27115, 16840, 27205, 7682,
                                                                       7742, 17460, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28195, 0, 3,
                                                                       27205, 16900, 27295, 7742,
                                                                       7802, 17560, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28345, 0, 3,
                                                                       27295, 16960, 27385, 7802,
                                                                       7862, 17660, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28495, 0, 3,
                                                                       27385, 17020, 27475, 7862,
                                                                       7922, 17760, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28645, 0, 3,
                                                                       27475, 17080, 27565, 7922,
                                                                       7982, 17860, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28795, 0, 3,
                                                                       27565, 17140, 27655, 7982,
                                                                       8042, 17960, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28945, 0, 3,
                                                                       27745, 17260, 27895, 8162,
                                                                       8252, 18060, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29170, 0, 3,
                                                                       27895, 17360, 28045, 8252,
                                                                       8342, 18210, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29395, 0, 3,
                                                                       28045, 17460, 28195, 8342,
                                                                       8432, 18360, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29620, 0, 3,
                                                                       28195, 17560, 28345, 8432,
                                                                       8522, 18510, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29845, 0, 3,
                                                                       28345, 17660, 28495, 8522,
                                                                       8612, 18660, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 30070, 0, 3,
                                                                       28495, 17760, 28645, 8612,
                                                                       8702, 18810, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 30295, 0, 3,
                                                                       28645, 17860, 28795, 8702,
                                                                       8792, 18960, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 30520, 0, 3,
                                                                       28945, 18060, 29170, 8972,
                                                                       9098, 19110, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 30835, 0, 3,
                                                                       29170, 18210, 29395, 9098,
                                                                       9224, 19320, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 31150, 0, 3,
                                                                       29395, 18360, 29620, 9224,
                                                                       9350, 19530, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 31465, 0, 3,
                                                                       29620, 18510, 29845, 9350,
                                                                       9476, 19740, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 31780, 0, 3,
                                                                       29845, 18660, 30070, 9476,
                                                                       9602, 19950, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 32095, 0, 3,
                                                                       30070, 18810, 30295, 9602,
                                                                       9728, 20160, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 32410, 0, 3,
                                                                       30520, 19110, 30835, 9980,
                                                                       10148, 20370, ncols,
                                                                       gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 32830, 0, 3,
                                                                       30835, 19320, 31150,
                                                                       10148, 10316, 20650,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 33250, 0, 3,
                                                                       31150, 19530, 31465,
                                                                       10316, 10484, 20930,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 33670, 0, 3,
                                                                       31465, 19740, 31780,
                                                                       10484, 10652, 21210,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 34090, 0, 3,
                                                                       31780, 19950, 32095,
                                                                       10652, 10820, 21490,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 34510, 0, 3,
                                                                       32410, 20370, 32830,
                                                                       11156, 11372, 21770,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 35050, 0, 3,
                                                                       32830, 20650, 33250,
                                                                       11372, 11588, 22130,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 35590, 0, 3,
                                                                       33250, 20930, 33670,
                                                                       11588, 11804, 22490,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 36130, 0, 3,
                                                                       33670, 21210, 34090,
                                                                       11804, 12020, 22850,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 36670, 0, 3,
                                                                       34510, 21770, 35050,
                                                                       12452, 12722, 23210,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 37345, 0, 3,
                                                                       35050, 22130, 35590,
                                                                       12722, 12992, 23660,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 38020, 0, 3,
                                                                       35590, 22490, 36130,
                                                                       12992, 13262, 24110,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 38695, 0, 3,
                                                                       36670, 23210, 37345,
                                                                       13802, 14132, 24560,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 39520, 0, 3,
                                                                       37345, 23660, 38020,
                                                                       14132, 14462, 25110,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 40345, 0, 3,
                                                                       38695, 24560, 39520,
                                                                       15122, 15518, 25660,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 41335, 30520, 315, ncols);

                    simdfunc::contract_primitives(buffer, 41839, 32410, 420, ncols);

                    simdfunc::contract_primitives(buffer, 42511, 34510, 540, ncols);

                    simdfunc::contract_primitives(buffer, 43375, 36670, 675, ncols);

                    simdfunc::contract_primitives(buffer, 44455, 38695, 825, ncols);

                    simdfunc::contract_primitives(buffer, 45775, 40345, 990, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 41650, 41335, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 42259, 41839, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 43051, 42511, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 44050, 43375, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 45280, 44455, 55, 1, nmax);

        simdtrf::transform_g_inner(buffer, 46765, 45775, 66, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 47359, 41650, 42259, 9, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 47926, 42259, 43051, 9, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 48682, 43051, 44050, 9, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 49654, 44050, 45280, 9, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 50869, 45280, 46765, 9, nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 52354, 47359, 47926, 9, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 53488, 47926, 48682, 9, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 55000, 48682, 49654, 9, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 56944, 49654, 50869, 9, nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 59374, 52354, 53488, 9, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 61264, 53488, 55000, 9, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 63784, 55000, 56944, 9, nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 67024, 59374, 61264, 9, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 69859, 61264, 63784, 9, nmax);

        simdtrf::compute_hrr_hh_out_of_first(buffer, coordinates, 73639, 67024, 69859, 9, nmax);

        simdtrf::transform_h_inner(buffer, 77608, 73639, 21, 9, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 77608, 99, nmax);
    }

    for (size_t m = 0; m < 1089; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
