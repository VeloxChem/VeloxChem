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


#include "SimdThreeCenterElectronRepulsionRecIHH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
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
compute_ihh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ihh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 154596, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1573 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 154596, 93920, 8998, dimensions);

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
                                                        16}, ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 24, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 69, 0, 3, 8, 9,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 75, 0, 3, 9, 10,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 81, 0, 3, 10, 11,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 87, 0, 3, 11, 12,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 93, 0, 3, 12, 13,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 99, 0, 3, 13, 14,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 105, 0, 3, 14, 15,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 111, 0, 3, 15, 16,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 117, 0, 3, 16, 17,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 123, 0, 3, 17, 18,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 129, 0, 3, 18, 19,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 135, 0, 3, 19, 20,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 141, 0, 3, 20, 21,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 147, 0, 3, 21, 22,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 24, 27,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 27, 30,
                                                                       75, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 30, 33,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 33, 36,
                                                                       87, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 36, 39,
                                                                       93, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 39, 42,
                                                                       99, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 42, 45,
                                                                       105, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 45, 48,
                                                                       111, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 48, 51,
                                                                       117, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 51, 54,
                                                                       123, 129, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 54, 57,
                                                                       129, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 263, 0, 3, 57, 60,
                                                                       135, 141, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 60, 63,
                                                                       141, 147, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 283, 0, 3, 69, 75,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 75, 81,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 313, 0, 3, 81, 87,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 87, 93,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 343, 0, 3, 93, 99,
                                                                       193, 203, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 99,
                                                                       105, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 373, 0, 3, 105,
                                                                       111, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 111,
                                                                       117, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 403, 0, 3, 117,
                                                                       123, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 123,
                                                                       129, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 433, 0, 3, 129,
                                                                       135, 253, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 135,
                                                                       141, 263, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 463, 0, 3, 153,
                                                                       163, 283, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 484, 0, 3, 163,
                                                                       173, 298, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 505, 0, 3, 173,
                                                                       183, 313, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 526, 0, 3, 183,
                                                                       193, 328, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 547, 0, 3, 193,
                                                                       203, 343, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 568, 0, 3, 203,
                                                                       213, 358, 373, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 589, 0, 3, 213,
                                                                       223, 373, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 610, 0, 3, 223,
                                                                       233, 388, 403, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 631, 0, 3, 233,
                                                                       243, 403, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 652, 0, 3, 243,
                                                                       253, 418, 433, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 673, 0, 3, 253,
                                                                       263, 433, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 694, 0, 3, 283,
                                                                       298, 463, 484, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 722, 0, 3, 298,
                                                                       313, 484, 505, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 750, 0, 3, 313,
                                                                       328, 505, 526, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 778, 0, 3, 328,
                                                                       343, 526, 547, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 806, 0, 3, 343,
                                                                       358, 547, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 834, 0, 3, 358,
                                                                       373, 568, 589, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 862, 0, 3, 373,
                                                                       388, 589, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 890, 0, 3, 388,
                                                                       403, 610, 631, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 918, 0, 3, 403,
                                                                       418, 631, 652, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 946, 0, 3, 418,
                                                                       433, 652, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 974, 0, 3, 463,
                                                                       484, 694, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1010, 0, 3, 484,
                                                                       505, 722, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1046, 0, 3, 505,
                                                                       526, 750, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1082, 0, 3, 526,
                                                                       547, 778, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1118, 0, 3, 547,
                                                                       568, 806, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1154, 0, 3, 568,
                                                                       589, 834, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1190, 0, 3, 589,
                                                                       610, 862, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1226, 0, 3, 610,
                                                                       631, 890, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1262, 0, 3, 631,
                                                                       652, 918, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1298, 0, 3, 694,
                                                                       722, 974, 1010, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1343, 0, 3, 722,
                                                                       750, 1010, 1046, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1388, 0, 3, 750,
                                                                       778, 1046, 1082, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1433, 0, 3, 778,
                                                                       806, 1082, 1118, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1478, 0, 3, 806,
                                                                       834, 1118, 1154, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1523, 0, 3, 834,
                                                                       862, 1154, 1190, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1568, 0, 3, 862,
                                                                       890, 1190, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1613, 0, 3, 890,
                                                                       918, 1226, 1262, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1658, 0, 3, 974,
                                                                       1010, 1298, 1343, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1713, 0, 3, 1010,
                                                                       1046, 1343, 1388, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1768, 0, 3, 1046,
                                                                       1082, 1388, 1433, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1823, 0, 3, 1082,
                                                                       1118, 1433, 1478, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1878, 0, 3, 1118,
                                                                       1154, 1478, 1523, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1933, 0, 3, 1154,
                                                                       1190, 1523, 1568, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1190,
                                                                       1226, 1568, 1613, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2043, 0, 3, 1298,
                                                                       1343, 1658, 1713, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2109, 0, 3, 1343,
                                                                       1388, 1713, 1768, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2175, 0, 3, 1388,
                                                                       1433, 1768, 1823, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2241, 0, 3, 1433,
                                                                       1478, 1823, 1878, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2307, 0, 3, 1478,
                                                                       1523, 1878, 1933, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2373, 0, 3, 1523,
                                                                       1568, 1933, 1988, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2439, 0, 3, 1658,
                                                                       1713, 2043, 2109, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2517, 0, 3, 1713,
                                                                       1768, 2109, 2175, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2595, 0, 3, 1768,
                                                                       1823, 2175, 2241, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2673, 0, 3, 1823,
                                                                       1878, 2241, 2307, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2751, 0, 3, 1878,
                                                                       1933, 2307, 2373, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2829, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2832, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2835, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2838, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2841, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2844, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2847, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2850, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2853, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2856, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2859, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2862, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2865, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2868, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2871, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2874, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2877, 3, 8, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2886, 3, 9, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2895, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2904, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2913, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2922, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2931, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2940, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2949, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2958, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2967, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2976, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2985, 3, 20, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2994, 3, 21, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3003, 3, 22, 66,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3012, 3, 24, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3030, 3, 27, 75,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3048, 3, 30, 81,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3066, 3, 33, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3084, 3, 36, 93,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3102, 3, 39, 99,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3120, 3, 42, 105,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3138, 3, 45, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3156, 3, 48, 117,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3174, 3, 51, 123,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3192, 3, 54, 129,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3210, 3, 57, 135,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3228, 3, 60, 141,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3246, 3, 63, 147,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3264, 3, 69, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3294, 3, 75, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3324, 3, 81, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3354, 3, 87, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3384, 3, 93, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3414, 3, 99, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3444, 3, 105, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3474, 3, 111, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3504, 3, 117, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3534, 3, 123, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3564, 3, 129, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3594, 3, 135, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3624, 3, 141, 273,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3654, 3, 153, 283,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3699, 3, 163, 298,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3744, 3, 173, 313,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3789, 3, 183, 328,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3834, 3, 193, 343,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3879, 3, 203, 358,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3924, 3, 213, 373,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3969, 3, 223, 388,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4014, 3, 233, 403,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4059, 3, 243, 418,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4104, 3, 253, 433,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4149, 3, 263, 448,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4194, 3, 283, 463,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4257, 3, 298, 484,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4320, 3, 313, 505,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4383, 3, 328, 526,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4446, 3, 343, 547,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4509, 3, 358, 568,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4572, 3, 373, 589,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4635, 3, 388, 610,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4698, 3, 403, 631,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4761, 3, 418, 652,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4824, 3, 433, 673,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4887, 3, 463, 694,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4971, 3, 484, 722,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5055, 3, 505, 750,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5139, 3, 526, 778,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5223, 3, 547, 806,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5307, 3, 568, 834,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5391, 3, 589, 862,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5475, 3, 610, 890,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5559, 3, 631, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5643, 3, 652, 946,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5727, 3, 694, 974,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5835, 3, 722,
                                                                       1010, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5943, 3, 750,
                                                                       1046, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6051, 3, 778,
                                                                       1082, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6159, 3, 806,
                                                                       1118, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6267, 3, 834,
                                                                       1154, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6375, 3, 862,
                                                                       1190, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6483, 3, 890,
                                                                       1226, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6591, 3, 918,
                                                                       1262, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6699, 3, 974,
                                                                       1298, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6834, 3, 1010,
                                                                       1343, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6969, 3, 1046,
                                                                       1388, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7104, 3, 1082,
                                                                       1433, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7239, 3, 1118,
                                                                       1478, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7374, 3, 1154,
                                                                       1523, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7509, 3, 1190,
                                                                       1568, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7644, 3, 1226,
                                                                       1613, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7779, 3, 1298,
                                                                       1658, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7944, 3, 1343,
                                                                       1713, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8109, 3, 1388,
                                                                       1768, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8274, 3, 1433,
                                                                       1823, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8439, 3, 1478,
                                                                       1878, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8604, 3, 1523,
                                                                       1933, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8769, 3, 1568,
                                                                       1988, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8934, 3, 1658,
                                                                       2043, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9132, 3, 1713,
                                                                       2109, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9330, 3, 1768,
                                                                       2175, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9528, 3, 1823,
                                                                       2241, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9726, 3, 1878,
                                                                       2307, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9924, 3, 1933,
                                                                       2373, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10122, 3, 2043,
                                                                       2439, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10356, 3, 2109,
                                                                       2517, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10590, 3, 2175,
                                                                       2595, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10824, 3, 2241,
                                                                       2673, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 11058, 3, 2307,
                                                                       2751, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11292, 3, 8, 9,
                                                                       2835, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11298, 3, 9, 10,
                                                                       2838, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11304, 3, 10, 11,
                                                                       2841, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11310, 3, 11, 12,
                                                                       2844, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11316, 3, 12, 13,
                                                                       2847, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11322, 3, 13, 14,
                                                                       2850, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11328, 3, 14, 15,
                                                                       2853, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11334, 3, 15, 16,
                                                                       2856, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11340, 3, 16, 17,
                                                                       2859, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11346, 3, 17, 18,
                                                                       2862, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11352, 3, 18, 19,
                                                                       2865, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11358, 3, 19, 20,
                                                                       2868, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11364, 3, 20, 21,
                                                                       2871, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11370, 3, 21, 22,
                                                                       2874, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11376, 0, 3,
                                                                       11292, 2835, 11298, 2895,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11394, 0, 3,
                                                                       11298, 2838, 11304, 2904,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11412, 0, 3,
                                                                       11304, 2841, 11310, 2913,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11430, 0, 3,
                                                                       11310, 2844, 11316, 2922,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11448, 0, 3,
                                                                       11316, 2847, 11322, 2931,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11466, 0, 3,
                                                                       11322, 2850, 11328, 2940,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11484, 0, 3,
                                                                       11328, 2853, 11334, 2949,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11502, 0, 3,
                                                                       11334, 2856, 11340, 2958,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11520, 0, 3,
                                                                       11340, 2859, 11346, 2967,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11538, 0, 3,
                                                                       11346, 2862, 11352, 2976,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11556, 0, 3,
                                                                       11352, 2865, 11358, 2985,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11574, 0, 3,
                                                                       11358, 2868, 11364, 2994,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11592, 0, 3,
                                                                       11364, 2871, 11370, 3003,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11610, 0, 3,
                                                                       11376, 2895, 11394, 69,
                                                                       75, 3048, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11646, 0, 3,
                                                                       11394, 2904, 11412, 75,
                                                                       81, 3066, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11682, 0, 3,
                                                                       11412, 2913, 11430, 81,
                                                                       87, 3084, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11718, 0, 3,
                                                                       11430, 2922, 11448, 87,
                                                                       93, 3102, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11754, 0, 3,
                                                                       11448, 2931, 11466, 93,
                                                                       99, 3120, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11790, 0, 3,
                                                                       11466, 2940, 11484, 99,
                                                                       105, 3138, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11826, 0, 3,
                                                                       11484, 2949, 11502, 105,
                                                                       111, 3156, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11862, 0, 3,
                                                                       11502, 2958, 11520, 111,
                                                                       117, 3174, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11898, 0, 3,
                                                                       11520, 2967, 11538, 117,
                                                                       123, 3192, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11934, 0, 3,
                                                                       11538, 2976, 11556, 123,
                                                                       129, 3210, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11970, 0, 3,
                                                                       11556, 2985, 11574, 129,
                                                                       135, 3228, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12006, 0, 3,
                                                                       11574, 2994, 11592, 135,
                                                                       141, 3246, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12042, 0, 3,
                                                                       11610, 3048, 11646, 153,
                                                                       163, 3324, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12102, 0, 3,
                                                                       11646, 3066, 11682, 163,
                                                                       173, 3354, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12162, 0, 3,
                                                                       11682, 3084, 11718, 173,
                                                                       183, 3384, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12222, 0, 3,
                                                                       11718, 3102, 11754, 183,
                                                                       193, 3414, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12282, 0, 3,
                                                                       11754, 3120, 11790, 193,
                                                                       203, 3444, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12342, 0, 3,
                                                                       11790, 3138, 11826, 203,
                                                                       213, 3474, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12402, 0, 3,
                                                                       11826, 3156, 11862, 213,
                                                                       223, 3504, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12462, 0, 3,
                                                                       11862, 3174, 11898, 223,
                                                                       233, 3534, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12522, 0, 3,
                                                                       11898, 3192, 11934, 233,
                                                                       243, 3564, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12582, 0, 3,
                                                                       11934, 3210, 11970, 243,
                                                                       253, 3594, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12642, 0, 3,
                                                                       11970, 3228, 12006, 253,
                                                                       263, 3624, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12702, 0, 3,
                                                                       12042, 3324, 12102, 283,
                                                                       298, 3744, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12792, 0, 3,
                                                                       12102, 3354, 12162, 298,
                                                                       313, 3789, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12882, 0, 3,
                                                                       12162, 3384, 12222, 313,
                                                                       328, 3834, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12972, 0, 3,
                                                                       12222, 3414, 12282, 328,
                                                                       343, 3879, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13062, 0, 3,
                                                                       12282, 3444, 12342, 343,
                                                                       358, 3924, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13152, 0, 3,
                                                                       12342, 3474, 12402, 358,
                                                                       373, 3969, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13242, 0, 3,
                                                                       12402, 3504, 12462, 373,
                                                                       388, 4014, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13332, 0, 3,
                                                                       12462, 3534, 12522, 388,
                                                                       403, 4059, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13422, 0, 3,
                                                                       12522, 3564, 12582, 403,
                                                                       418, 4104, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13512, 0, 3,
                                                                       12582, 3594, 12642, 418,
                                                                       433, 4149, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13602, 0, 3,
                                                                       12702, 3744, 12792, 463,
                                                                       484, 4320, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13728, 0, 3,
                                                                       12792, 3789, 12882, 484,
                                                                       505, 4383, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13854, 0, 3,
                                                                       12882, 3834, 12972, 505,
                                                                       526, 4446, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13980, 0, 3,
                                                                       12972, 3879, 13062, 526,
                                                                       547, 4509, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14106, 0, 3,
                                                                       13062, 3924, 13152, 547,
                                                                       568, 4572, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14232, 0, 3,
                                                                       13152, 3969, 13242, 568,
                                                                       589, 4635, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14358, 0, 3,
                                                                       13242, 4014, 13332, 589,
                                                                       610, 4698, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14484, 0, 3,
                                                                       13332, 4059, 13422, 610,
                                                                       631, 4761, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14610, 0, 3,
                                                                       13422, 4104, 13512, 631,
                                                                       652, 4824, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14736, 0, 3,
                                                                       13602, 4320, 13728, 694,
                                                                       722, 5055, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14904, 0, 3,
                                                                       13728, 4383, 13854, 722,
                                                                       750, 5139, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15072, 0, 3,
                                                                       13854, 4446, 13980, 750,
                                                                       778, 5223, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15240, 0, 3,
                                                                       13980, 4509, 14106, 778,
                                                                       806, 5307, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15408, 0, 3,
                                                                       14106, 4572, 14232, 806,
                                                                       834, 5391, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15576, 0, 3,
                                                                       14232, 4635, 14358, 834,
                                                                       862, 5475, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15744, 0, 3,
                                                                       14358, 4698, 14484, 862,
                                                                       890, 5559, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15912, 0, 3,
                                                                       14484, 4761, 14610, 890,
                                                                       918, 5643, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16080, 0, 3,
                                                                       14736, 5055, 14904, 974,
                                                                       1010, 5943, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16296, 0, 3,
                                                                       14904, 5139, 15072, 1010,
                                                                       1046, 6051, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16512, 0, 3,
                                                                       15072, 5223, 15240, 1046,
                                                                       1082, 6159, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16728, 0, 3,
                                                                       15240, 5307, 15408, 1082,
                                                                       1118, 6267, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16944, 0, 3,
                                                                       15408, 5391, 15576, 1118,
                                                                       1154, 6375, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17160, 0, 3,
                                                                       15576, 5475, 15744, 1154,
                                                                       1190, 6483, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17376, 0, 3,
                                                                       15744, 5559, 15912, 1190,
                                                                       1226, 6591, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17592, 0, 3,
                                                                       16080, 5943, 16296, 1298,
                                                                       1343, 6969, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17862, 0, 3,
                                                                       16296, 6051, 16512, 1343,
                                                                       1388, 7104, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18132, 0, 3,
                                                                       16512, 6159, 16728, 1388,
                                                                       1433, 7239, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18402, 0, 3,
                                                                       16728, 6267, 16944, 1433,
                                                                       1478, 7374, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18672, 0, 3,
                                                                       16944, 6375, 17160, 1478,
                                                                       1523, 7509, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18942, 0, 3,
                                                                       17160, 6483, 17376, 1523,
                                                                       1568, 7644, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19212, 0, 3,
                                                                       17592, 6969, 17862, 1658,
                                                                       1713, 8109, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19542, 0, 3,
                                                                       17862, 7104, 18132, 1713,
                                                                       1768, 8274, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19872, 0, 3,
                                                                       18132, 7239, 18402, 1768,
                                                                       1823, 8439, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20202, 0, 3,
                                                                       18402, 7374, 18672, 1823,
                                                                       1878, 8604, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20532, 0, 3,
                                                                       18672, 7509, 18942, 1878,
                                                                       1933, 8769, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 20862, 0, 3,
                                                                       19212, 8109, 19542, 2043,
                                                                       2109, 9330, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 21258, 0, 3,
                                                                       19542, 8274, 19872, 2109,
                                                                       2175, 9528, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 21654, 0, 3,
                                                                       19872, 8439, 20202, 2175,
                                                                       2241, 9726, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 22050, 0, 3,
                                                                       20202, 8604, 20532, 2241,
                                                                       2307, 9924, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 22446, 0, 3,
                                                                       20862, 9330, 21258, 2439,
                                                                       2517, 10590, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 22914, 0, 3,
                                                                       21258, 9528, 21654, 2517,
                                                                       2595, 10824, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 23382, 0, 3,
                                                                       21654, 9726, 22050, 2595,
                                                                       2673, 11058, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23850, 3, 2829,
                                                                       2832, 11292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23860, 3, 2832,
                                                                       2835, 11298, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23870, 3, 2835,
                                                                       2838, 11304, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23880, 3, 2838,
                                                                       2841, 11310, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23890, 3, 2841,
                                                                       2844, 11316, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23900, 3, 2844,
                                                                       2847, 11322, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23910, 3, 2847,
                                                                       2850, 11328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23920, 3, 2850,
                                                                       2853, 11334, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23930, 3, 2853,
                                                                       2856, 11340, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23940, 3, 2856,
                                                                       2859, 11346, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23950, 3, 2859,
                                                                       2862, 11352, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23960, 3, 2862,
                                                                       2865, 11358, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23970, 3, 2865,
                                                                       2868, 11364, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23980, 3, 2868,
                                                                       2871, 11370, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 23990, 0, 3,
                                                                       23850, 11292, 23860, 2877,
                                                                       2886, 11376, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24020, 0, 3,
                                                                       23860, 11298, 23870, 2886,
                                                                       2895, 11394, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24050, 0, 3,
                                                                       23870, 11304, 23880, 2895,
                                                                       2904, 11412, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24080, 0, 3,
                                                                       23880, 11310, 23890, 2904,
                                                                       2913, 11430, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24110, 0, 3,
                                                                       23890, 11316, 23900, 2913,
                                                                       2922, 11448, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24140, 0, 3,
                                                                       23900, 11322, 23910, 2922,
                                                                       2931, 11466, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24170, 0, 3,
                                                                       23910, 11328, 23920, 2931,
                                                                       2940, 11484, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24200, 0, 3,
                                                                       23920, 11334, 23930, 2940,
                                                                       2949, 11502, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24230, 0, 3,
                                                                       23930, 11340, 23940, 2949,
                                                                       2958, 11520, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24260, 0, 3,
                                                                       23940, 11346, 23950, 2958,
                                                                       2967, 11538, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24290, 0, 3,
                                                                       23950, 11352, 23960, 2967,
                                                                       2976, 11556, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24320, 0, 3,
                                                                       23960, 11358, 23970, 2976,
                                                                       2985, 11574, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24350, 0, 3,
                                                                       23970, 11364, 23980, 2985,
                                                                       2994, 11592, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24380, 0, 3,
                                                                       23990, 11376, 24020, 3012,
                                                                       3030, 11610, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24440, 0, 3,
                                                                       24020, 11394, 24050, 3030,
                                                                       3048, 11646, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24500, 0, 3,
                                                                       24050, 11412, 24080, 3048,
                                                                       3066, 11682, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24560, 0, 3,
                                                                       24080, 11430, 24110, 3066,
                                                                       3084, 11718, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24620, 0, 3,
                                                                       24110, 11448, 24140, 3084,
                                                                       3102, 11754, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24680, 0, 3,
                                                                       24140, 11466, 24170, 3102,
                                                                       3120, 11790, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24740, 0, 3,
                                                                       24170, 11484, 24200, 3120,
                                                                       3138, 11826, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24800, 0, 3,
                                                                       24200, 11502, 24230, 3138,
                                                                       3156, 11862, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24860, 0, 3,
                                                                       24230, 11520, 24260, 3156,
                                                                       3174, 11898, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24920, 0, 3,
                                                                       24260, 11538, 24290, 3174,
                                                                       3192, 11934, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24980, 0, 3,
                                                                       24290, 11556, 24320, 3192,
                                                                       3210, 11970, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25040, 0, 3,
                                                                       24320, 11574, 24350, 3210,
                                                                       3228, 12006, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25100, 0, 3,
                                                                       24380, 11610, 24440, 3264,
                                                                       3294, 12042, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25200, 0, 3,
                                                                       24440, 11646, 24500, 3294,
                                                                       3324, 12102, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25300, 0, 3,
                                                                       24500, 11682, 24560, 3324,
                                                                       3354, 12162, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25400, 0, 3,
                                                                       24560, 11718, 24620, 3354,
                                                                       3384, 12222, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25500, 0, 3,
                                                                       24620, 11754, 24680, 3384,
                                                                       3414, 12282, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25600, 0, 3,
                                                                       24680, 11790, 24740, 3414,
                                                                       3444, 12342, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25700, 0, 3,
                                                                       24740, 11826, 24800, 3444,
                                                                       3474, 12402, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25800, 0, 3,
                                                                       24800, 11862, 24860, 3474,
                                                                       3504, 12462, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25900, 0, 3,
                                                                       24860, 11898, 24920, 3504,
                                                                       3534, 12522, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26000, 0, 3,
                                                                       24920, 11934, 24980, 3534,
                                                                       3564, 12582, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26100, 0, 3,
                                                                       24980, 11970, 25040, 3564,
                                                                       3594, 12642, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26200, 0, 3,
                                                                       25100, 12042, 25200, 3654,
                                                                       3699, 12702, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26350, 0, 3,
                                                                       25200, 12102, 25300, 3699,
                                                                       3744, 12792, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26500, 0, 3,
                                                                       25300, 12162, 25400, 3744,
                                                                       3789, 12882, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26650, 0, 3,
                                                                       25400, 12222, 25500, 3789,
                                                                       3834, 12972, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26800, 0, 3,
                                                                       25500, 12282, 25600, 3834,
                                                                       3879, 13062, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26950, 0, 3,
                                                                       25600, 12342, 25700, 3879,
                                                                       3924, 13152, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27100, 0, 3,
                                                                       25700, 12402, 25800, 3924,
                                                                       3969, 13242, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27250, 0, 3,
                                                                       25800, 12462, 25900, 3969,
                                                                       4014, 13332, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27400, 0, 3,
                                                                       25900, 12522, 26000, 4014,
                                                                       4059, 13422, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27550, 0, 3,
                                                                       26000, 12582, 26100, 4059,
                                                                       4104, 13512, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 27700, 0, 3,
                                                                       26200, 12702, 26350, 4194,
                                                                       4257, 13602, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 27910, 0, 3,
                                                                       26350, 12792, 26500, 4257,
                                                                       4320, 13728, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28120, 0, 3,
                                                                       26500, 12882, 26650, 4320,
                                                                       4383, 13854, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28330, 0, 3,
                                                                       26650, 12972, 26800, 4383,
                                                                       4446, 13980, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28540, 0, 3,
                                                                       26800, 13062, 26950, 4446,
                                                                       4509, 14106, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28750, 0, 3,
                                                                       26950, 13152, 27100, 4509,
                                                                       4572, 14232, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28960, 0, 3,
                                                                       27100, 13242, 27250, 4572,
                                                                       4635, 14358, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29170, 0, 3,
                                                                       27250, 13332, 27400, 4635,
                                                                       4698, 14484, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29380, 0, 3,
                                                                       27400, 13422, 27550, 4698,
                                                                       4761, 14610, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 29590, 0, 3,
                                                                       27700, 13602, 27910, 4887,
                                                                       4971, 14736, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 29870, 0, 3,
                                                                       27910, 13728, 28120, 4971,
                                                                       5055, 14904, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30150, 0, 3,
                                                                       28120, 13854, 28330, 5055,
                                                                       5139, 15072, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30430, 0, 3,
                                                                       28330, 13980, 28540, 5139,
                                                                       5223, 15240, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30710, 0, 3,
                                                                       28540, 14106, 28750, 5223,
                                                                       5307, 15408, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30990, 0, 3,
                                                                       28750, 14232, 28960, 5307,
                                                                       5391, 15576, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 31270, 0, 3,
                                                                       28960, 14358, 29170, 5391,
                                                                       5475, 15744, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 31550, 0, 3,
                                                                       29170, 14484, 29380, 5475,
                                                                       5559, 15912, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 31830, 0, 3,
                                                                       29590, 14736, 29870, 5727,
                                                                       5835, 16080, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32190, 0, 3,
                                                                       29870, 14904, 30150, 5835,
                                                                       5943, 16296, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32550, 0, 3,
                                                                       30150, 15072, 30430, 5943,
                                                                       6051, 16512, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32910, 0, 3,
                                                                       30430, 15240, 30710, 6051,
                                                                       6159, 16728, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33270, 0, 3,
                                                                       30710, 15408, 30990, 6159,
                                                                       6267, 16944, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33630, 0, 3,
                                                                       30990, 15576, 31270, 6267,
                                                                       6375, 17160, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33990, 0, 3,
                                                                       31270, 15744, 31550, 6375,
                                                                       6483, 17376, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 34350, 0, 3,
                                                                       31830, 16080, 32190, 6699,
                                                                       6834, 17592, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 34800, 0, 3,
                                                                       32190, 16296, 32550, 6834,
                                                                       6969, 17862, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 35250, 0, 3,
                                                                       32550, 16512, 32910, 6969,
                                                                       7104, 18132, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 35700, 0, 3,
                                                                       32910, 16728, 33270, 7104,
                                                                       7239, 18402, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 36150, 0, 3,
                                                                       33270, 16944, 33630, 7239,
                                                                       7374, 18672, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 36600, 0, 3,
                                                                       33630, 17160, 33990, 7374,
                                                                       7509, 18942, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 37050, 0, 3,
                                                                       34350, 17592, 34800, 7779,
                                                                       7944, 19212, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 37600, 0, 3,
                                                                       34800, 17862, 35250, 7944,
                                                                       8109, 19542, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 38150, 0, 3,
                                                                       35250, 18132, 35700, 8109,
                                                                       8274, 19872, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 38700, 0, 3,
                                                                       35700, 18402, 36150, 8274,
                                                                       8439, 20202, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 39250, 0, 3,
                                                                       36150, 18672, 36600, 8439,
                                                                       8604, 20532, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 39800, 0, 3,
                                                                       37050, 19212, 37600, 8934,
                                                                       9132, 20862, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 40460, 0, 3,
                                                                       37600, 19542, 38150, 9132,
                                                                       9330, 21258, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 41120, 0, 3,
                                                                       38150, 19872, 38700, 9330,
                                                                       9528, 21654, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 41780, 0, 3,
                                                                       38700, 20202, 39250, 9528,
                                                                       9726, 22050, ncols, gamma,
                                                                       p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 42440, 0, 3,
                                                                       39800, 20862, 40460,
                                                                       10122, 10356, 22446,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 43220, 0, 3,
                                                                       40460, 21258, 41120,
                                                                       10356, 10590, 22914,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 44000, 0, 3,
                                                                       41120, 21654, 41780,
                                                                       10590, 10824, 23382,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44780, 3, 11292,
                                                                       11298, 23870, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44795, 3, 11298,
                                                                       11304, 23880, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44810, 3, 11304,
                                                                       11310, 23890, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44825, 3, 11310,
                                                                       11316, 23900, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44840, 3, 11316,
                                                                       11322, 23910, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44855, 3, 11322,
                                                                       11328, 23920, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44870, 3, 11328,
                                                                       11334, 23930, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44885, 3, 11334,
                                                                       11340, 23940, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44900, 3, 11340,
                                                                       11346, 23950, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44915, 3, 11346,
                                                                       11352, 23960, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44930, 3, 11352,
                                                                       11358, 23970, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44945, 3, 11358,
                                                                       11364, 23980, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 44960, 0, 3,
                                                                       44780, 23870, 44795,
                                                                       11376, 11394, 24050,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45005, 0, 3,
                                                                       44795, 23880, 44810,
                                                                       11394, 11412, 24080,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45050, 0, 3,
                                                                       44810, 23890, 44825,
                                                                       11412, 11430, 24110,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45095, 0, 3,
                                                                       44825, 23900, 44840,
                                                                       11430, 11448, 24140,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45140, 0, 3,
                                                                       44840, 23910, 44855,
                                                                       11448, 11466, 24170,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45185, 0, 3,
                                                                       44855, 23920, 44870,
                                                                       11466, 11484, 24200,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45230, 0, 3,
                                                                       44870, 23930, 44885,
                                                                       11484, 11502, 24230,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45275, 0, 3,
                                                                       44885, 23940, 44900,
                                                                       11502, 11520, 24260,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45320, 0, 3,
                                                                       44900, 23950, 44915,
                                                                       11520, 11538, 24290,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45365, 0, 3,
                                                                       44915, 23960, 44930,
                                                                       11538, 11556, 24320,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45410, 0, 3,
                                                                       44930, 23970, 44945,
                                                                       11556, 11574, 24350,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45455, 0, 3,
                                                                       44960, 24050, 45005,
                                                                       11610, 11646, 24500,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45545, 0, 3,
                                                                       45005, 24080, 45050,
                                                                       11646, 11682, 24560,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45635, 0, 3,
                                                                       45050, 24110, 45095,
                                                                       11682, 11718, 24620,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45725, 0, 3,
                                                                       45095, 24140, 45140,
                                                                       11718, 11754, 24680,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45815, 0, 3,
                                                                       45140, 24170, 45185,
                                                                       11754, 11790, 24740,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45905, 0, 3,
                                                                       45185, 24200, 45230,
                                                                       11790, 11826, 24800,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45995, 0, 3,
                                                                       45230, 24230, 45275,
                                                                       11826, 11862, 24860,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46085, 0, 3,
                                                                       45275, 24260, 45320,
                                                                       11862, 11898, 24920,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46175, 0, 3,
                                                                       45320, 24290, 45365,
                                                                       11898, 11934, 24980,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46265, 0, 3,
                                                                       45365, 24320, 45410,
                                                                       11934, 11970, 25040,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46355, 0, 3,
                                                                       45455, 24500, 45545,
                                                                       12042, 12102, 25300,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46505, 0, 3,
                                                                       45545, 24560, 45635,
                                                                       12102, 12162, 25400,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46655, 0, 3,
                                                                       45635, 24620, 45725,
                                                                       12162, 12222, 25500,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46805, 0, 3,
                                                                       45725, 24680, 45815,
                                                                       12222, 12282, 25600,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46955, 0, 3,
                                                                       45815, 24740, 45905,
                                                                       12282, 12342, 25700,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47105, 0, 3,
                                                                       45905, 24800, 45995,
                                                                       12342, 12402, 25800,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47255, 0, 3,
                                                                       45995, 24860, 46085,
                                                                       12402, 12462, 25900,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47405, 0, 3,
                                                                       46085, 24920, 46175,
                                                                       12462, 12522, 26000,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47555, 0, 3,
                                                                       46175, 24980, 46265,
                                                                       12522, 12582, 26100,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 47705, 0, 3,
                                                                       46355, 25300, 46505,
                                                                       12702, 12792, 26500,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 47930, 0, 3,
                                                                       46505, 25400, 46655,
                                                                       12792, 12882, 26650,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48155, 0, 3,
                                                                       46655, 25500, 46805,
                                                                       12882, 12972, 26800,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48380, 0, 3,
                                                                       46805, 25600, 46955,
                                                                       12972, 13062, 26950,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48605, 0, 3,
                                                                       46955, 25700, 47105,
                                                                       13062, 13152, 27100,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48830, 0, 3,
                                                                       47105, 25800, 47255,
                                                                       13152, 13242, 27250,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49055, 0, 3,
                                                                       47255, 25900, 47405,
                                                                       13242, 13332, 27400,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49280, 0, 3,
                                                                       47405, 26000, 47555,
                                                                       13332, 13422, 27550,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 49505, 0, 3,
                                                                       47705, 26500, 47930,
                                                                       13602, 13728, 28120,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 49820, 0, 3,
                                                                       47930, 26650, 48155,
                                                                       13728, 13854, 28330,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50135, 0, 3,
                                                                       48155, 26800, 48380,
                                                                       13854, 13980, 28540,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50450, 0, 3,
                                                                       48380, 26950, 48605,
                                                                       13980, 14106, 28750,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50765, 0, 3,
                                                                       48605, 27100, 48830,
                                                                       14106, 14232, 28960,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51080, 0, 3,
                                                                       48830, 27250, 49055,
                                                                       14232, 14358, 29170,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51395, 0, 3,
                                                                       49055, 27400, 49280,
                                                                       14358, 14484, 29380,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 51710, 0, 3,
                                                                       49505, 28120, 49820,
                                                                       14736, 14904, 30150,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 52130, 0, 3,
                                                                       49820, 28330, 50135,
                                                                       14904, 15072, 30430,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 52550, 0, 3,
                                                                       50135, 28540, 50450,
                                                                       15072, 15240, 30710,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 52970, 0, 3,
                                                                       50450, 28750, 50765,
                                                                       15240, 15408, 30990,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 53390, 0, 3,
                                                                       50765, 28960, 51080,
                                                                       15408, 15576, 31270,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 53810, 0, 3,
                                                                       51080, 29170, 51395,
                                                                       15576, 15744, 31550,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 54230, 0, 3,
                                                                       51710, 30150, 52130,
                                                                       16080, 16296, 32550,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 54770, 0, 3,
                                                                       52130, 30430, 52550,
                                                                       16296, 16512, 32910,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 55310, 0, 3,
                                                                       52550, 30710, 52970,
                                                                       16512, 16728, 33270,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 55850, 0, 3,
                                                                       52970, 30990, 53390,
                                                                       16728, 16944, 33630,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 56390, 0, 3,
                                                                       53390, 31270, 53810,
                                                                       16944, 17160, 33990,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 56930, 0, 3,
                                                                       54230, 32550, 54770,
                                                                       17592, 17862, 35250,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 57605, 0, 3,
                                                                       54770, 32910, 55310,
                                                                       17862, 18132, 35700,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 58280, 0, 3,
                                                                       55310, 33270, 55850,
                                                                       18132, 18402, 36150,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 58955, 0, 3,
                                                                       55850, 33630, 56390,
                                                                       18402, 18672, 36600,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 59630, 0, 3,
                                                                       56930, 35250, 57605,
                                                                       19212, 19542, 38150,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 60455, 0, 3,
                                                                       57605, 35700, 58280,
                                                                       19542, 19872, 38700,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 61280, 0, 3,
                                                                       58280, 36150, 58955,
                                                                       19872, 20202, 39250,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 62105, 0, 3,
                                                                       59630, 38150, 60455,
                                                                       20862, 21258, 41120,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 63095, 0, 3,
                                                                       60455, 38700, 61280,
                                                                       21258, 21654, 41780,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 64085, 0, 3,
                                                                       62105, 41120, 63095,
                                                                       22446, 22914, 44000,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65255, 3, 23850,
                                                                       23860, 44780, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65276, 3, 23860,
                                                                       23870, 44795, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65297, 3, 23870,
                                                                       23880, 44810, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65318, 3, 23880,
                                                                       23890, 44825, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65339, 3, 23890,
                                                                       23900, 44840, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65360, 3, 23900,
                                                                       23910, 44855, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65381, 3, 23910,
                                                                       23920, 44870, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65402, 3, 23920,
                                                                       23930, 44885, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65423, 3, 23930,
                                                                       23940, 44900, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65444, 3, 23940,
                                                                       23950, 44915, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65465, 3, 23950,
                                                                       23960, 44930, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65486, 3, 23960,
                                                                       23970, 44945, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65507, 0, 3,
                                                                       65255, 44780, 65276,
                                                                       23990, 24020, 44960,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65570, 0, 3,
                                                                       65276, 44795, 65297,
                                                                       24020, 24050, 45005,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65633, 0, 3,
                                                                       65297, 44810, 65318,
                                                                       24050, 24080, 45050,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65696, 0, 3,
                                                                       65318, 44825, 65339,
                                                                       24080, 24110, 45095,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65759, 0, 3,
                                                                       65339, 44840, 65360,
                                                                       24110, 24140, 45140,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65822, 0, 3,
                                                                       65360, 44855, 65381,
                                                                       24140, 24170, 45185,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65885, 0, 3,
                                                                       65381, 44870, 65402,
                                                                       24170, 24200, 45230,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65948, 0, 3,
                                                                       65402, 44885, 65423,
                                                                       24200, 24230, 45275,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66011, 0, 3,
                                                                       65423, 44900, 65444,
                                                                       24230, 24260, 45320,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66074, 0, 3,
                                                                       65444, 44915, 65465,
                                                                       24260, 24290, 45365,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66137, 0, 3,
                                                                       65465, 44930, 65486,
                                                                       24290, 24320, 45410,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66200, 0, 3,
                                                                       65507, 44960, 65570,
                                                                       24380, 24440, 45455,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66326, 0, 3,
                                                                       65570, 45005, 65633,
                                                                       24440, 24500, 45545,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66452, 0, 3,
                                                                       65633, 45050, 65696,
                                                                       24500, 24560, 45635,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66578, 0, 3,
                                                                       65696, 45095, 65759,
                                                                       24560, 24620, 45725,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66704, 0, 3,
                                                                       65759, 45140, 65822,
                                                                       24620, 24680, 45815,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66830, 0, 3,
                                                                       65822, 45185, 65885,
                                                                       24680, 24740, 45905,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66956, 0, 3,
                                                                       65885, 45230, 65948,
                                                                       24740, 24800, 45995,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67082, 0, 3,
                                                                       65948, 45275, 66011,
                                                                       24800, 24860, 46085,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67208, 0, 3,
                                                                       66011, 45320, 66074,
                                                                       24860, 24920, 46175,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67334, 0, 3,
                                                                       66074, 45365, 66137,
                                                                       24920, 24980, 46265,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 67460, 0, 3,
                                                                       66200, 45455, 66326,
                                                                       25100, 25200, 46355,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 67670, 0, 3,
                                                                       66326, 45545, 66452,
                                                                       25200, 25300, 46505,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 67880, 0, 3,
                                                                       66452, 45635, 66578,
                                                                       25300, 25400, 46655,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68090, 0, 3,
                                                                       66578, 45725, 66704,
                                                                       25400, 25500, 46805,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68300, 0, 3,
                                                                       66704, 45815, 66830,
                                                                       25500, 25600, 46955,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68510, 0, 3,
                                                                       66830, 45905, 66956,
                                                                       25600, 25700, 47105,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68720, 0, 3,
                                                                       66956, 45995, 67082,
                                                                       25700, 25800, 47255,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68930, 0, 3,
                                                                       67082, 46085, 67208,
                                                                       25800, 25900, 47405,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 69140, 0, 3,
                                                                       67208, 46175, 67334,
                                                                       25900, 26000, 47555,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 69350, 0, 3,
                                                                       67460, 46355, 67670,
                                                                       26200, 26350, 47705,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 69665, 0, 3,
                                                                       67670, 46505, 67880,
                                                                       26350, 26500, 47930,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 69980, 0, 3,
                                                                       67880, 46655, 68090,
                                                                       26500, 26650, 48155,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 70295, 0, 3,
                                                                       68090, 46805, 68300,
                                                                       26650, 26800, 48380,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 70610, 0, 3,
                                                                       68300, 46955, 68510,
                                                                       26800, 26950, 48605,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 70925, 0, 3,
                                                                       68510, 47105, 68720,
                                                                       26950, 27100, 48830,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 71240, 0, 3,
                                                                       68720, 47255, 68930,
                                                                       27100, 27250, 49055,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 71555, 0, 3,
                                                                       68930, 47405, 69140,
                                                                       27250, 27400, 49280,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 71870, 0, 3,
                                                                       69350, 47705, 69665,
                                                                       27700, 27910, 49505,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 72311, 0, 3,
                                                                       69665, 47930, 69980,
                                                                       27910, 28120, 49820,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 72752, 0, 3,
                                                                       69980, 48155, 70295,
                                                                       28120, 28330, 50135,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 73193, 0, 3,
                                                                       70295, 48380, 70610,
                                                                       28330, 28540, 50450,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 73634, 0, 3,
                                                                       70610, 48605, 70925,
                                                                       28540, 28750, 50765,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 74075, 0, 3,
                                                                       70925, 48830, 71240,
                                                                       28750, 28960, 51080,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 74516, 0, 3,
                                                                       71240, 49055, 71555,
                                                                       28960, 29170, 51395,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 74957, 0, 3,
                                                                       71870, 49505, 72311,
                                                                       29590, 29870, 51710,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 75545, 0, 3,
                                                                       72311, 49820, 72752,
                                                                       29870, 30150, 52130,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 76133, 0, 3,
                                                                       72752, 50135, 73193,
                                                                       30150, 30430, 52550,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 76721, 0, 3,
                                                                       73193, 50450, 73634,
                                                                       30430, 30710, 52970,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 77309, 0, 3,
                                                                       73634, 50765, 74075,
                                                                       30710, 30990, 53390,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 77897, 0, 3,
                                                                       74075, 51080, 74516,
                                                                       30990, 31270, 53810,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 78485, 0, 3,
                                                                       74957, 51710, 75545,
                                                                       31830, 32190, 54230,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 79241, 0, 3,
                                                                       75545, 52130, 76133,
                                                                       32190, 32550, 54770,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 79997, 0, 3,
                                                                       76133, 52550, 76721,
                                                                       32550, 32910, 55310,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 80753, 0, 3,
                                                                       76721, 52970, 77309,
                                                                       32910, 33270, 55850,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 81509, 0, 3,
                                                                       77309, 53390, 77897,
                                                                       33270, 33630, 56390,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 82265, 0, 3,
                                                                       78485, 54230, 79241,
                                                                       34350, 34800, 56930,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 83210, 0, 3,
                                                                       79241, 54770, 79997,
                                                                       34800, 35250, 57605,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 84155, 0, 3,
                                                                       79997, 55310, 80753,
                                                                       35250, 35700, 58280,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 85100, 0, 3,
                                                                       80753, 55850, 81509,
                                                                       35700, 36150, 58955,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 86045, 0, 3,
                                                                       82265, 56930, 83210,
                                                                       37050, 37600, 59630,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 87200, 0, 3,
                                                                       83210, 57605, 84155,
                                                                       37600, 38150, 60455,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 88355, 0, 3,
                                                                       84155, 58280, 85100,
                                                                       38150, 38700, 61280,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 89510, 0, 3,
                                                                       86045, 59630, 87200,
                                                                       39800, 40460, 62105,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 90896, 0, 3,
                                                                       87200, 60455, 88355,
                                                                       40460, 41120, 63095,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 92282, 0, 3,
                                                                       89510, 62105, 90896,
                                                                       42440, 43220, 64085,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 93920, 74957, 588, ncols);

                    simdfunc::contract_primitives(buffer, 94816, 78485, 756, ncols);

                    simdfunc::contract_primitives(buffer, 95968, 82265, 945, ncols);

                    simdfunc::contract_primitives(buffer, 97408, 86045, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 99168, 89510, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 101280, 92282, 1638, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 94508, 93920, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 95572, 94816, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 96913, 95968, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 98563, 97408, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 100554, 99168, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 102918, 101280, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 103776, 94508, 95572, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 104700, 95572, 96913, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 105888, 96913, 98563, 11,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 107373, 98563, 100554, 11,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 109188, 100554, 102918, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 111366, 103776, 104700, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 113214, 104700, 105888, 11,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 115590, 105888, 107373, 11,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 118560, 107373, 109188, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 122190, 111366, 113214, 11,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 125270, 113214, 115590, 11,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 129230, 115590, 118560, 11,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 134180, 122190, 125270, 11,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 138800, 125270, 129230, 11,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 144740, 134180, 138800, 11,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 151208, 144740, 28, 11, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 151208, 121, nmax);
    }

    for (size_t m = 0; m < 1573; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
