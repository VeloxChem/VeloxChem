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


#include "SimdThreeCenterElectronRepulsionRecHGK.hpp"

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
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hgk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hgk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 158853, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1485 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 158853, 121008, 8610, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2043, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2046, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2049, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2052, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2055, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2058, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2061, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2064, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2067, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2070, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2073, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2076, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2079, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2082, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2085, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2088, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2091, 3, 8, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2100, 3, 9, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2109, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2118, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2127, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2136, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2145, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2154, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2163, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2172, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2181, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2190, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2199, 3, 20, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2208, 3, 21, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2217, 3, 22, 66,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2226, 3, 24, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2244, 3, 27, 75,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2262, 3, 30, 81,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2280, 3, 33, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2298, 3, 36, 93,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2316, 3, 39, 99,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2334, 3, 42, 105,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2352, 3, 45, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2370, 3, 48, 117,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2388, 3, 51, 123,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2406, 3, 54, 129,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2424, 3, 57, 135,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2442, 3, 60, 141,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2460, 3, 63, 147,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2478, 3, 69, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2508, 3, 75, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2538, 3, 81, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2568, 3, 87, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2598, 3, 93, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2628, 3, 99, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2658, 3, 105, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2688, 3, 111, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2718, 3, 117, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2748, 3, 123, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2778, 3, 129, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2808, 3, 135, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2838, 3, 141, 273,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2868, 3, 153, 283,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2913, 3, 163, 298,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2958, 3, 173, 313,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3003, 3, 183, 328,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3048, 3, 193, 343,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3093, 3, 203, 358,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3138, 3, 213, 373,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3183, 3, 223, 388,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3228, 3, 233, 403,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3273, 3, 243, 418,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3318, 3, 253, 433,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3363, 3, 263, 448,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3408, 3, 283, 463,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3471, 3, 298, 484,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3534, 3, 313, 505,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3597, 3, 328, 526,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3660, 3, 343, 547,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3723, 3, 358, 568,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3786, 3, 373, 589,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3849, 3, 388, 610,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3912, 3, 403, 631,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3975, 3, 418, 652,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4038, 3, 433, 673,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4101, 3, 463, 694,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4185, 3, 484, 722,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4269, 3, 505, 750,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4353, 3, 526, 778,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4437, 3, 547, 806,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4521, 3, 568, 834,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4605, 3, 589, 862,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4689, 3, 610, 890,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4773, 3, 631, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4857, 3, 652, 946,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4941, 3, 694, 974,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5049, 3, 722,
                                                                       1010, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5157, 3, 750,
                                                                       1046, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5265, 3, 778,
                                                                       1082, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5373, 3, 806,
                                                                       1118, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5481, 3, 834,
                                                                       1154, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5589, 3, 862,
                                                                       1190, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5697, 3, 890,
                                                                       1226, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5805, 3, 918,
                                                                       1262, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5913, 3, 974,
                                                                       1298, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6048, 3, 1010,
                                                                       1343, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6183, 3, 1046,
                                                                       1388, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6318, 3, 1082,
                                                                       1433, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6453, 3, 1118,
                                                                       1478, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6588, 3, 1154,
                                                                       1523, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6723, 3, 1190,
                                                                       1568, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6858, 3, 1226,
                                                                       1613, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6993, 3, 1298,
                                                                       1658, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7158, 3, 1343,
                                                                       1713, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7323, 3, 1388,
                                                                       1768, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7488, 3, 1433,
                                                                       1823, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7653, 3, 1478,
                                                                       1878, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7818, 3, 1523,
                                                                       1933, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7983, 3, 1568,
                                                                       1988, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8148, 3, 8, 9,
                                                                       2049, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8154, 3, 9, 10,
                                                                       2052, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8160, 3, 10, 11,
                                                                       2055, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8166, 3, 11, 12,
                                                                       2058, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8172, 3, 12, 13,
                                                                       2061, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8178, 3, 13, 14,
                                                                       2064, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8184, 3, 14, 15,
                                                                       2067, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8190, 3, 15, 16,
                                                                       2070, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8196, 3, 16, 17,
                                                                       2073, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8202, 3, 17, 18,
                                                                       2076, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8208, 3, 18, 19,
                                                                       2079, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8214, 3, 19, 20,
                                                                       2082, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8220, 3, 20, 21,
                                                                       2085, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8226, 3, 21, 22,
                                                                       2088, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8232, 0, 3, 8148,
                                                                       2049, 8154, 2109, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8250, 0, 3, 8154,
                                                                       2052, 8160, 2118, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8268, 0, 3, 8160,
                                                                       2055, 8166, 2127, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8286, 0, 3, 8166,
                                                                       2058, 8172, 2136, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8304, 0, 3, 8172,
                                                                       2061, 8178, 2145, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8322, 0, 3, 8178,
                                                                       2064, 8184, 2154, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8340, 0, 3, 8184,
                                                                       2067, 8190, 2163, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8358, 0, 3, 8190,
                                                                       2070, 8196, 2172, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8376, 0, 3, 8196,
                                                                       2073, 8202, 2181, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8394, 0, 3, 8202,
                                                                       2076, 8208, 2190, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8412, 0, 3, 8208,
                                                                       2079, 8214, 2199, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8430, 0, 3, 8214,
                                                                       2082, 8220, 2208, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8448, 0, 3, 8220,
                                                                       2085, 8226, 2217, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8466, 0, 3, 8232,
                                                                       2109, 8250, 69, 75, 2262,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8502, 0, 3, 8250,
                                                                       2118, 8268, 75, 81, 2280,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8538, 0, 3, 8268,
                                                                       2127, 8286, 81, 87, 2298,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8574, 0, 3, 8286,
                                                                       2136, 8304, 87, 93, 2316,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8610, 0, 3, 8304,
                                                                       2145, 8322, 93, 99, 2334,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8646, 0, 3, 8322,
                                                                       2154, 8340, 99, 105, 2352,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8682, 0, 3, 8340,
                                                                       2163, 8358, 105, 111,
                                                                       2370, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8718, 0, 3, 8358,
                                                                       2172, 8376, 111, 117,
                                                                       2388, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8754, 0, 3, 8376,
                                                                       2181, 8394, 117, 123,
                                                                       2406, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8790, 0, 3, 8394,
                                                                       2190, 8412, 123, 129,
                                                                       2424, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8826, 0, 3, 8412,
                                                                       2199, 8430, 129, 135,
                                                                       2442, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8862, 0, 3, 8430,
                                                                       2208, 8448, 135, 141,
                                                                       2460, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8898, 0, 3, 8466,
                                                                       2262, 8502, 153, 163,
                                                                       2538, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8958, 0, 3, 8502,
                                                                       2280, 8538, 163, 173,
                                                                       2568, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9018, 0, 3, 8538,
                                                                       2298, 8574, 173, 183,
                                                                       2598, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9078, 0, 3, 8574,
                                                                       2316, 8610, 183, 193,
                                                                       2628, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9138, 0, 3, 8610,
                                                                       2334, 8646, 193, 203,
                                                                       2658, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9198, 0, 3, 8646,
                                                                       2352, 8682, 203, 213,
                                                                       2688, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9258, 0, 3, 8682,
                                                                       2370, 8718, 213, 223,
                                                                       2718, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9318, 0, 3, 8718,
                                                                       2388, 8754, 223, 233,
                                                                       2748, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9378, 0, 3, 8754,
                                                                       2406, 8790, 233, 243,
                                                                       2778, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9438, 0, 3, 8790,
                                                                       2424, 8826, 243, 253,
                                                                       2808, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9498, 0, 3, 8826,
                                                                       2442, 8862, 253, 263,
                                                                       2838, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9558, 0, 3, 8898,
                                                                       2538, 8958, 283, 298,
                                                                       2958, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9648, 0, 3, 8958,
                                                                       2568, 9018, 298, 313,
                                                                       3003, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9738, 0, 3, 9018,
                                                                       2598, 9078, 313, 328,
                                                                       3048, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9828, 0, 3, 9078,
                                                                       2628, 9138, 328, 343,
                                                                       3093, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9918, 0, 3, 9138,
                                                                       2658, 9198, 343, 358,
                                                                       3138, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10008, 0, 3, 9198,
                                                                       2688, 9258, 358, 373,
                                                                       3183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10098, 0, 3, 9258,
                                                                       2718, 9318, 373, 388,
                                                                       3228, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10188, 0, 3, 9318,
                                                                       2748, 9378, 388, 403,
                                                                       3273, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10278, 0, 3, 9378,
                                                                       2778, 9438, 403, 418,
                                                                       3318, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10368, 0, 3, 9438,
                                                                       2808, 9498, 418, 433,
                                                                       3363, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10458, 0, 3, 9558,
                                                                       2958, 9648, 463, 484,
                                                                       3534, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10584, 0, 3, 9648,
                                                                       3003, 9738, 484, 505,
                                                                       3597, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10710, 0, 3, 9738,
                                                                       3048, 9828, 505, 526,
                                                                       3660, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10836, 0, 3, 9828,
                                                                       3093, 9918, 526, 547,
                                                                       3723, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10962, 0, 3, 9918,
                                                                       3138, 10008, 547, 568,
                                                                       3786, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11088, 0, 3,
                                                                       10008, 3183, 10098, 568,
                                                                       589, 3849, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11214, 0, 3,
                                                                       10098, 3228, 10188, 589,
                                                                       610, 3912, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11340, 0, 3,
                                                                       10188, 3273, 10278, 610,
                                                                       631, 3975, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11466, 0, 3,
                                                                       10278, 3318, 10368, 631,
                                                                       652, 4038, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11592, 0, 3,
                                                                       10458, 3534, 10584, 694,
                                                                       722, 4269, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11760, 0, 3,
                                                                       10584, 3597, 10710, 722,
                                                                       750, 4353, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11928, 0, 3,
                                                                       10710, 3660, 10836, 750,
                                                                       778, 4437, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12096, 0, 3,
                                                                       10836, 3723, 10962, 778,
                                                                       806, 4521, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12264, 0, 3,
                                                                       10962, 3786, 11088, 806,
                                                                       834, 4605, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12432, 0, 3,
                                                                       11088, 3849, 11214, 834,
                                                                       862, 4689, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12600, 0, 3,
                                                                       11214, 3912, 11340, 862,
                                                                       890, 4773, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12768, 0, 3,
                                                                       11340, 3975, 11466, 890,
                                                                       918, 4857, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12936, 0, 3,
                                                                       11592, 4269, 11760, 974,
                                                                       1010, 5157, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13152, 0, 3,
                                                                       11760, 4353, 11928, 1010,
                                                                       1046, 5265, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13368, 0, 3,
                                                                       11928, 4437, 12096, 1046,
                                                                       1082, 5373, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13584, 0, 3,
                                                                       12096, 4521, 12264, 1082,
                                                                       1118, 5481, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13800, 0, 3,
                                                                       12264, 4605, 12432, 1118,
                                                                       1154, 5589, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14016, 0, 3,
                                                                       12432, 4689, 12600, 1154,
                                                                       1190, 5697, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14232, 0, 3,
                                                                       12600, 4773, 12768, 1190,
                                                                       1226, 5805, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14448, 0, 3,
                                                                       12936, 5157, 13152, 1298,
                                                                       1343, 6183, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14718, 0, 3,
                                                                       13152, 5265, 13368, 1343,
                                                                       1388, 6318, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14988, 0, 3,
                                                                       13368, 5373, 13584, 1388,
                                                                       1433, 6453, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 15258, 0, 3,
                                                                       13584, 5481, 13800, 1433,
                                                                       1478, 6588, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 15528, 0, 3,
                                                                       13800, 5589, 14016, 1478,
                                                                       1523, 6723, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 15798, 0, 3,
                                                                       14016, 5697, 14232, 1523,
                                                                       1568, 6858, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 16068, 0, 3,
                                                                       14448, 6183, 14718, 1658,
                                                                       1713, 7323, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 16398, 0, 3,
                                                                       14718, 6318, 14988, 1713,
                                                                       1768, 7488, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 16728, 0, 3,
                                                                       14988, 6453, 15258, 1768,
                                                                       1823, 7653, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 17058, 0, 3,
                                                                       15258, 6588, 15528, 1823,
                                                                       1878, 7818, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 17388, 0, 3,
                                                                       15528, 6723, 15798, 1878,
                                                                       1933, 7983, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17718, 3, 2043,
                                                                       2046, 8148, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17728, 3, 2046,
                                                                       2049, 8154, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17738, 3, 2049,
                                                                       2052, 8160, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17748, 3, 2052,
                                                                       2055, 8166, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17758, 3, 2055,
                                                                       2058, 8172, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17768, 3, 2058,
                                                                       2061, 8178, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17778, 3, 2061,
                                                                       2064, 8184, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17788, 3, 2064,
                                                                       2067, 8190, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17798, 3, 2067,
                                                                       2070, 8196, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17808, 3, 2070,
                                                                       2073, 8202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17818, 3, 2073,
                                                                       2076, 8208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17828, 3, 2076,
                                                                       2079, 8214, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17838, 3, 2079,
                                                                       2082, 8220, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17848, 3, 2082,
                                                                       2085, 8226, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17858, 0, 3,
                                                                       17718, 8148, 17728, 2091,
                                                                       2100, 8232, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17888, 0, 3,
                                                                       17728, 8154, 17738, 2100,
                                                                       2109, 8250, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17918, 0, 3,
                                                                       17738, 8160, 17748, 2109,
                                                                       2118, 8268, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17948, 0, 3,
                                                                       17748, 8166, 17758, 2118,
                                                                       2127, 8286, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17978, 0, 3,
                                                                       17758, 8172, 17768, 2127,
                                                                       2136, 8304, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18008, 0, 3,
                                                                       17768, 8178, 17778, 2136,
                                                                       2145, 8322, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18038, 0, 3,
                                                                       17778, 8184, 17788, 2145,
                                                                       2154, 8340, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18068, 0, 3,
                                                                       17788, 8190, 17798, 2154,
                                                                       2163, 8358, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18098, 0, 3,
                                                                       17798, 8196, 17808, 2163,
                                                                       2172, 8376, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18128, 0, 3,
                                                                       17808, 8202, 17818, 2172,
                                                                       2181, 8394, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18158, 0, 3,
                                                                       17818, 8208, 17828, 2181,
                                                                       2190, 8412, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18188, 0, 3,
                                                                       17828, 8214, 17838, 2190,
                                                                       2199, 8430, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18218, 0, 3,
                                                                       17838, 8220, 17848, 2199,
                                                                       2208, 8448, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18248, 0, 3,
                                                                       17858, 8232, 17888, 2226,
                                                                       2244, 8466, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18308, 0, 3,
                                                                       17888, 8250, 17918, 2244,
                                                                       2262, 8502, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18368, 0, 3,
                                                                       17918, 8268, 17948, 2262,
                                                                       2280, 8538, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18428, 0, 3,
                                                                       17948, 8286, 17978, 2280,
                                                                       2298, 8574, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18488, 0, 3,
                                                                       17978, 8304, 18008, 2298,
                                                                       2316, 8610, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18548, 0, 3,
                                                                       18008, 8322, 18038, 2316,
                                                                       2334, 8646, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18608, 0, 3,
                                                                       18038, 8340, 18068, 2334,
                                                                       2352, 8682, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18668, 0, 3,
                                                                       18068, 8358, 18098, 2352,
                                                                       2370, 8718, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18728, 0, 3,
                                                                       18098, 8376, 18128, 2370,
                                                                       2388, 8754, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18788, 0, 3,
                                                                       18128, 8394, 18158, 2388,
                                                                       2406, 8790, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18848, 0, 3,
                                                                       18158, 8412, 18188, 2406,
                                                                       2424, 8826, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18908, 0, 3,
                                                                       18188, 8430, 18218, 2424,
                                                                       2442, 8862, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18968, 0, 3,
                                                                       18248, 8466, 18308, 2478,
                                                                       2508, 8898, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19068, 0, 3,
                                                                       18308, 8502, 18368, 2508,
                                                                       2538, 8958, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19168, 0, 3,
                                                                       18368, 8538, 18428, 2538,
                                                                       2568, 9018, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19268, 0, 3,
                                                                       18428, 8574, 18488, 2568,
                                                                       2598, 9078, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19368, 0, 3,
                                                                       18488, 8610, 18548, 2598,
                                                                       2628, 9138, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19468, 0, 3,
                                                                       18548, 8646, 18608, 2628,
                                                                       2658, 9198, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19568, 0, 3,
                                                                       18608, 8682, 18668, 2658,
                                                                       2688, 9258, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19668, 0, 3,
                                                                       18668, 8718, 18728, 2688,
                                                                       2718, 9318, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19768, 0, 3,
                                                                       18728, 8754, 18788, 2718,
                                                                       2748, 9378, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19868, 0, 3,
                                                                       18788, 8790, 18848, 2748,
                                                                       2778, 9438, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19968, 0, 3,
                                                                       18848, 8826, 18908, 2778,
                                                                       2808, 9498, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20068, 0, 3,
                                                                       18968, 8898, 19068, 2868,
                                                                       2913, 9558, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20218, 0, 3,
                                                                       19068, 8958, 19168, 2913,
                                                                       2958, 9648, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20368, 0, 3,
                                                                       19168, 9018, 19268, 2958,
                                                                       3003, 9738, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20518, 0, 3,
                                                                       19268, 9078, 19368, 3003,
                                                                       3048, 9828, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20668, 0, 3,
                                                                       19368, 9138, 19468, 3048,
                                                                       3093, 9918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20818, 0, 3,
                                                                       19468, 9198, 19568, 3093,
                                                                       3138, 10008, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20968, 0, 3,
                                                                       19568, 9258, 19668, 3138,
                                                                       3183, 10098, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21118, 0, 3,
                                                                       19668, 9318, 19768, 3183,
                                                                       3228, 10188, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21268, 0, 3,
                                                                       19768, 9378, 19868, 3228,
                                                                       3273, 10278, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21418, 0, 3,
                                                                       19868, 9438, 19968, 3273,
                                                                       3318, 10368, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 21568, 0, 3,
                                                                       20068, 9558, 20218, 3408,
                                                                       3471, 10458, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 21778, 0, 3,
                                                                       20218, 9648, 20368, 3471,
                                                                       3534, 10584, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 21988, 0, 3,
                                                                       20368, 9738, 20518, 3534,
                                                                       3597, 10710, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22198, 0, 3,
                                                                       20518, 9828, 20668, 3597,
                                                                       3660, 10836, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22408, 0, 3,
                                                                       20668, 9918, 20818, 3660,
                                                                       3723, 10962, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22618, 0, 3,
                                                                       20818, 10008, 20968, 3723,
                                                                       3786, 11088, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22828, 0, 3,
                                                                       20968, 10098, 21118, 3786,
                                                                       3849, 11214, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23038, 0, 3,
                                                                       21118, 10188, 21268, 3849,
                                                                       3912, 11340, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23248, 0, 3,
                                                                       21268, 10278, 21418, 3912,
                                                                       3975, 11466, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 23458, 0, 3,
                                                                       21568, 10458, 21778, 4101,
                                                                       4185, 11592, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 23738, 0, 3,
                                                                       21778, 10584, 21988, 4185,
                                                                       4269, 11760, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24018, 0, 3,
                                                                       21988, 10710, 22198, 4269,
                                                                       4353, 11928, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24298, 0, 3,
                                                                       22198, 10836, 22408, 4353,
                                                                       4437, 12096, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24578, 0, 3,
                                                                       22408, 10962, 22618, 4437,
                                                                       4521, 12264, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24858, 0, 3,
                                                                       22618, 11088, 22828, 4521,
                                                                       4605, 12432, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 25138, 0, 3,
                                                                       22828, 11214, 23038, 4605,
                                                                       4689, 12600, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 25418, 0, 3,
                                                                       23038, 11340, 23248, 4689,
                                                                       4773, 12768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 25698, 0, 3,
                                                                       23458, 11592, 23738, 4941,
                                                                       5049, 12936, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 26058, 0, 3,
                                                                       23738, 11760, 24018, 5049,
                                                                       5157, 13152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 26418, 0, 3,
                                                                       24018, 11928, 24298, 5157,
                                                                       5265, 13368, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 26778, 0, 3,
                                                                       24298, 12096, 24578, 5265,
                                                                       5373, 13584, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27138, 0, 3,
                                                                       24578, 12264, 24858, 5373,
                                                                       5481, 13800, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27498, 0, 3,
                                                                       24858, 12432, 25138, 5481,
                                                                       5589, 14016, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27858, 0, 3,
                                                                       25138, 12600, 25418, 5589,
                                                                       5697, 14232, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 28218, 0, 3,
                                                                       25698, 12936, 26058, 5913,
                                                                       6048, 14448, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 28668, 0, 3,
                                                                       26058, 13152, 26418, 6048,
                                                                       6183, 14718, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 29118, 0, 3,
                                                                       26418, 13368, 26778, 6183,
                                                                       6318, 14988, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 29568, 0, 3,
                                                                       26778, 13584, 27138, 6318,
                                                                       6453, 15258, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 30018, 0, 3,
                                                                       27138, 13800, 27498, 6453,
                                                                       6588, 15528, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 30468, 0, 3,
                                                                       27498, 14016, 27858, 6588,
                                                                       6723, 15798, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 30918, 0, 3,
                                                                       28218, 14448, 28668, 6993,
                                                                       7158, 16068, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 31468, 0, 3,
                                                                       28668, 14718, 29118, 7158,
                                                                       7323, 16398, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 32018, 0, 3,
                                                                       29118, 14988, 29568, 7323,
                                                                       7488, 16728, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 32568, 0, 3,
                                                                       29568, 15258, 30018, 7488,
                                                                       7653, 17058, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 33118, 0, 3,
                                                                       30018, 15528, 30468, 7653,
                                                                       7818, 17388, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33668, 3, 8148,
                                                                       8154, 17738, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33683, 3, 8154,
                                                                       8160, 17748, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33698, 3, 8160,
                                                                       8166, 17758, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33713, 3, 8166,
                                                                       8172, 17768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33728, 3, 8172,
                                                                       8178, 17778, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33743, 3, 8178,
                                                                       8184, 17788, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33758, 3, 8184,
                                                                       8190, 17798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33773, 3, 8190,
                                                                       8196, 17808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33788, 3, 8196,
                                                                       8202, 17818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33803, 3, 8202,
                                                                       8208, 17828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33818, 3, 8208,
                                                                       8214, 17838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33833, 3, 8214,
                                                                       8220, 17848, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 33848, 0, 3,
                                                                       33668, 17738, 33683, 8232,
                                                                       8250, 17918, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 33893, 0, 3,
                                                                       33683, 17748, 33698, 8250,
                                                                       8268, 17948, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 33938, 0, 3,
                                                                       33698, 17758, 33713, 8268,
                                                                       8286, 17978, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 33983, 0, 3,
                                                                       33713, 17768, 33728, 8286,
                                                                       8304, 18008, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34028, 0, 3,
                                                                       33728, 17778, 33743, 8304,
                                                                       8322, 18038, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34073, 0, 3,
                                                                       33743, 17788, 33758, 8322,
                                                                       8340, 18068, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34118, 0, 3,
                                                                       33758, 17798, 33773, 8340,
                                                                       8358, 18098, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34163, 0, 3,
                                                                       33773, 17808, 33788, 8358,
                                                                       8376, 18128, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34208, 0, 3,
                                                                       33788, 17818, 33803, 8376,
                                                                       8394, 18158, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34253, 0, 3,
                                                                       33803, 17828, 33818, 8394,
                                                                       8412, 18188, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34298, 0, 3,
                                                                       33818, 17838, 33833, 8412,
                                                                       8430, 18218, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34343, 0, 3,
                                                                       33848, 17918, 33893, 8466,
                                                                       8502, 18368, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34433, 0, 3,
                                                                       33893, 17948, 33938, 8502,
                                                                       8538, 18428, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34523, 0, 3,
                                                                       33938, 17978, 33983, 8538,
                                                                       8574, 18488, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34613, 0, 3,
                                                                       33983, 18008, 34028, 8574,
                                                                       8610, 18548, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34703, 0, 3,
                                                                       34028, 18038, 34073, 8610,
                                                                       8646, 18608, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34793, 0, 3,
                                                                       34073, 18068, 34118, 8646,
                                                                       8682, 18668, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34883, 0, 3,
                                                                       34118, 18098, 34163, 8682,
                                                                       8718, 18728, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34973, 0, 3,
                                                                       34163, 18128, 34208, 8718,
                                                                       8754, 18788, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35063, 0, 3,
                                                                       34208, 18158, 34253, 8754,
                                                                       8790, 18848, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35153, 0, 3,
                                                                       34253, 18188, 34298, 8790,
                                                                       8826, 18908, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35243, 0, 3,
                                                                       34343, 18368, 34433, 8898,
                                                                       8958, 19168, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35393, 0, 3,
                                                                       34433, 18428, 34523, 8958,
                                                                       9018, 19268, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35543, 0, 3,
                                                                       34523, 18488, 34613, 9018,
                                                                       9078, 19368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35693, 0, 3,
                                                                       34613, 18548, 34703, 9078,
                                                                       9138, 19468, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35843, 0, 3,
                                                                       34703, 18608, 34793, 9138,
                                                                       9198, 19568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35993, 0, 3,
                                                                       34793, 18668, 34883, 9198,
                                                                       9258, 19668, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 36143, 0, 3,
                                                                       34883, 18728, 34973, 9258,
                                                                       9318, 19768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 36293, 0, 3,
                                                                       34973, 18788, 35063, 9318,
                                                                       9378, 19868, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 36443, 0, 3,
                                                                       35063, 18848, 35153, 9378,
                                                                       9438, 19968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 36593, 0, 3,
                                                                       35243, 19168, 35393, 9558,
                                                                       9648, 20368, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 36818, 0, 3,
                                                                       35393, 19268, 35543, 9648,
                                                                       9738, 20518, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37043, 0, 3,
                                                                       35543, 19368, 35693, 9738,
                                                                       9828, 20668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37268, 0, 3,
                                                                       35693, 19468, 35843, 9828,
                                                                       9918, 20818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37493, 0, 3,
                                                                       35843, 19568, 35993, 9918,
                                                                       10008, 20968, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37718, 0, 3,
                                                                       35993, 19668, 36143,
                                                                       10008, 10098, 21118,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37943, 0, 3,
                                                                       36143, 19768, 36293,
                                                                       10098, 10188, 21268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 38168, 0, 3,
                                                                       36293, 19868, 36443,
                                                                       10188, 10278, 21418,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 38393, 0, 3,
                                                                       36593, 20368, 36818,
                                                                       10458, 10584, 21988,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 38708, 0, 3,
                                                                       36818, 20518, 37043,
                                                                       10584, 10710, 22198,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 39023, 0, 3,
                                                                       37043, 20668, 37268,
                                                                       10710, 10836, 22408,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 39338, 0, 3,
                                                                       37268, 20818, 37493,
                                                                       10836, 10962, 22618,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 39653, 0, 3,
                                                                       37493, 20968, 37718,
                                                                       10962, 11088, 22828,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 39968, 0, 3,
                                                                       37718, 21118, 37943,
                                                                       11088, 11214, 23038,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 40283, 0, 3,
                                                                       37943, 21268, 38168,
                                                                       11214, 11340, 23248,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 40598, 0, 3,
                                                                       38393, 21988, 38708,
                                                                       11592, 11760, 24018,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 41018, 0, 3,
                                                                       38708, 22198, 39023,
                                                                       11760, 11928, 24298,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 41438, 0, 3,
                                                                       39023, 22408, 39338,
                                                                       11928, 12096, 24578,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 41858, 0, 3,
                                                                       39338, 22618, 39653,
                                                                       12096, 12264, 24858,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 42278, 0, 3,
                                                                       39653, 22828, 39968,
                                                                       12264, 12432, 25138,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 42698, 0, 3,
                                                                       39968, 23038, 40283,
                                                                       12432, 12600, 25418,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 43118, 0, 3,
                                                                       40598, 24018, 41018,
                                                                       12936, 13152, 26418,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 43658, 0, 3,
                                                                       41018, 24298, 41438,
                                                                       13152, 13368, 26778,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 44198, 0, 3,
                                                                       41438, 24578, 41858,
                                                                       13368, 13584, 27138,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 44738, 0, 3,
                                                                       41858, 24858, 42278,
                                                                       13584, 13800, 27498,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 45278, 0, 3,
                                                                       42278, 25138, 42698,
                                                                       13800, 14016, 27858,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 45818, 0, 3,
                                                                       43118, 26418, 43658,
                                                                       14448, 14718, 29118,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 46493, 0, 3,
                                                                       43658, 26778, 44198,
                                                                       14718, 14988, 29568,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 47168, 0, 3,
                                                                       44198, 27138, 44738,
                                                                       14988, 15258, 30018,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 47843, 0, 3,
                                                                       44738, 27498, 45278,
                                                                       15258, 15528, 30468,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 48518, 0, 3,
                                                                       45818, 29118, 46493,
                                                                       16068, 16398, 32018,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 49343, 0, 3,
                                                                       46493, 29568, 47168,
                                                                       16398, 16728, 32568,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 50168, 0, 3,
                                                                       47168, 30018, 47843,
                                                                       16728, 17058, 33118,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50993, 3, 17718,
                                                                       17728, 33668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51014, 3, 17728,
                                                                       17738, 33683, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51035, 3, 17738,
                                                                       17748, 33698, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51056, 3, 17748,
                                                                       17758, 33713, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51077, 3, 17758,
                                                                       17768, 33728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51098, 3, 17768,
                                                                       17778, 33743, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51119, 3, 17778,
                                                                       17788, 33758, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51140, 3, 17788,
                                                                       17798, 33773, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51161, 3, 17798,
                                                                       17808, 33788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51182, 3, 17808,
                                                                       17818, 33803, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51203, 3, 17818,
                                                                       17828, 33818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51224, 3, 17828,
                                                                       17838, 33833, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51245, 0, 3,
                                                                       50993, 33668, 51014,
                                                                       17858, 17888, 33848,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51308, 0, 3,
                                                                       51014, 33683, 51035,
                                                                       17888, 17918, 33893,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51371, 0, 3,
                                                                       51035, 33698, 51056,
                                                                       17918, 17948, 33938,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51434, 0, 3,
                                                                       51056, 33713, 51077,
                                                                       17948, 17978, 33983,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51497, 0, 3,
                                                                       51077, 33728, 51098,
                                                                       17978, 18008, 34028,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51560, 0, 3,
                                                                       51098, 33743, 51119,
                                                                       18008, 18038, 34073,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51623, 0, 3,
                                                                       51119, 33758, 51140,
                                                                       18038, 18068, 34118,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51686, 0, 3,
                                                                       51140, 33773, 51161,
                                                                       18068, 18098, 34163,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51749, 0, 3,
                                                                       51161, 33788, 51182,
                                                                       18098, 18128, 34208,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51812, 0, 3,
                                                                       51182, 33803, 51203,
                                                                       18128, 18158, 34253,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51875, 0, 3,
                                                                       51203, 33818, 51224,
                                                                       18158, 18188, 34298,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 51938, 0, 3,
                                                                       51245, 33848, 51308,
                                                                       18248, 18308, 34343,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52064, 0, 3,
                                                                       51308, 33893, 51371,
                                                                       18308, 18368, 34433,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52190, 0, 3,
                                                                       51371, 33938, 51434,
                                                                       18368, 18428, 34523,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52316, 0, 3,
                                                                       51434, 33983, 51497,
                                                                       18428, 18488, 34613,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52442, 0, 3,
                                                                       51497, 34028, 51560,
                                                                       18488, 18548, 34703,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52568, 0, 3,
                                                                       51560, 34073, 51623,
                                                                       18548, 18608, 34793,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52694, 0, 3,
                                                                       51623, 34118, 51686,
                                                                       18608, 18668, 34883,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52820, 0, 3,
                                                                       51686, 34163, 51749,
                                                                       18668, 18728, 34973,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52946, 0, 3,
                                                                       51749, 34208, 51812,
                                                                       18728, 18788, 35063,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53072, 0, 3,
                                                                       51812, 34253, 51875,
                                                                       18788, 18848, 35153,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 53198, 0, 3,
                                                                       51938, 34343, 52064,
                                                                       18968, 19068, 35243,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 53408, 0, 3,
                                                                       52064, 34433, 52190,
                                                                       19068, 19168, 35393,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 53618, 0, 3,
                                                                       52190, 34523, 52316,
                                                                       19168, 19268, 35543,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 53828, 0, 3,
                                                                       52316, 34613, 52442,
                                                                       19268, 19368, 35693,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54038, 0, 3,
                                                                       52442, 34703, 52568,
                                                                       19368, 19468, 35843,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54248, 0, 3,
                                                                       52568, 34793, 52694,
                                                                       19468, 19568, 35993,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54458, 0, 3,
                                                                       52694, 34883, 52820,
                                                                       19568, 19668, 36143,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54668, 0, 3,
                                                                       52820, 34973, 52946,
                                                                       19668, 19768, 36293,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54878, 0, 3,
                                                                       52946, 35063, 53072,
                                                                       19768, 19868, 36443,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 55088, 0, 3,
                                                                       53198, 35243, 53408,
                                                                       20068, 20218, 36593,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 55403, 0, 3,
                                                                       53408, 35393, 53618,
                                                                       20218, 20368, 36818,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 55718, 0, 3,
                                                                       53618, 35543, 53828,
                                                                       20368, 20518, 37043,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 56033, 0, 3,
                                                                       53828, 35693, 54038,
                                                                       20518, 20668, 37268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 56348, 0, 3,
                                                                       54038, 35843, 54248,
                                                                       20668, 20818, 37493,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 56663, 0, 3,
                                                                       54248, 35993, 54458,
                                                                       20818, 20968, 37718,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 56978, 0, 3,
                                                                       54458, 36143, 54668,
                                                                       20968, 21118, 37943,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 57293, 0, 3,
                                                                       54668, 36293, 54878,
                                                                       21118, 21268, 38168,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 57608, 0, 3,
                                                                       55088, 36593, 55403,
                                                                       21568, 21778, 38393,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 58049, 0, 3,
                                                                       55403, 36818, 55718,
                                                                       21778, 21988, 38708,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 58490, 0, 3,
                                                                       55718, 37043, 56033,
                                                                       21988, 22198, 39023,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 58931, 0, 3,
                                                                       56033, 37268, 56348,
                                                                       22198, 22408, 39338,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 59372, 0, 3,
                                                                       56348, 37493, 56663,
                                                                       22408, 22618, 39653,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 59813, 0, 3,
                                                                       56663, 37718, 56978,
                                                                       22618, 22828, 39968,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 60254, 0, 3,
                                                                       56978, 37943, 57293,
                                                                       22828, 23038, 40283,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 60695, 0, 3,
                                                                       57608, 38393, 58049,
                                                                       23458, 23738, 40598,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 61283, 0, 3,
                                                                       58049, 38708, 58490,
                                                                       23738, 24018, 41018,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 61871, 0, 3,
                                                                       58490, 39023, 58931,
                                                                       24018, 24298, 41438,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 62459, 0, 3,
                                                                       58931, 39338, 59372,
                                                                       24298, 24578, 41858,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 63047, 0, 3,
                                                                       59372, 39653, 59813,
                                                                       24578, 24858, 42278,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 63635, 0, 3,
                                                                       59813, 39968, 60254,
                                                                       24858, 25138, 42698,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 64223, 0, 3,
                                                                       60695, 40598, 61283,
                                                                       25698, 26058, 43118,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 64979, 0, 3,
                                                                       61283, 41018, 61871,
                                                                       26058, 26418, 43658,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 65735, 0, 3,
                                                                       61871, 41438, 62459,
                                                                       26418, 26778, 44198,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 66491, 0, 3,
                                                                       62459, 41858, 63047,
                                                                       26778, 27138, 44738,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 67247, 0, 3,
                                                                       63047, 42278, 63635,
                                                                       27138, 27498, 45278,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 68003, 0, 3,
                                                                       64223, 43118, 64979,
                                                                       28218, 28668, 45818,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 68948, 0, 3,
                                                                       64979, 43658, 65735,
                                                                       28668, 29118, 46493,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 69893, 0, 3,
                                                                       65735, 44198, 66491,
                                                                       29118, 29568, 47168,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 70838, 0, 3,
                                                                       66491, 44738, 67247,
                                                                       29568, 30018, 47843,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 71783, 0, 3,
                                                                       68003, 45818, 68948,
                                                                       30918, 31468, 48518,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 72938, 0, 3,
                                                                       68948, 46493, 69893,
                                                                       31468, 32018, 49343,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 74093, 0, 3,
                                                                       69893, 47168, 70838,
                                                                       32018, 32568, 50168,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75248, 3, 33668,
                                                                       33683, 51035, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75276, 3, 33683,
                                                                       33698, 51056, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75304, 3, 33698,
                                                                       33713, 51077, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75332, 3, 33713,
                                                                       33728, 51098, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75360, 3, 33728,
                                                                       33743, 51119, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75388, 3, 33743,
                                                                       33758, 51140, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75416, 3, 33758,
                                                                       33773, 51161, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75444, 3, 33773,
                                                                       33788, 51182, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75472, 3, 33788,
                                                                       33803, 51203, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75500, 3, 33803,
                                                                       33818, 51224, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75528, 0, 3,
                                                                       75248, 51035, 75276,
                                                                       33848, 33893, 51371,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75612, 0, 3,
                                                                       75276, 51056, 75304,
                                                                       33893, 33938, 51434,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75696, 0, 3,
                                                                       75304, 51077, 75332,
                                                                       33938, 33983, 51497,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75780, 0, 3,
                                                                       75332, 51098, 75360,
                                                                       33983, 34028, 51560,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75864, 0, 3,
                                                                       75360, 51119, 75388,
                                                                       34028, 34073, 51623,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75948, 0, 3,
                                                                       75388, 51140, 75416,
                                                                       34073, 34118, 51686,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76032, 0, 3,
                                                                       75416, 51161, 75444,
                                                                       34118, 34163, 51749,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76116, 0, 3,
                                                                       75444, 51182, 75472,
                                                                       34163, 34208, 51812,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76200, 0, 3,
                                                                       75472, 51203, 75500,
                                                                       34208, 34253, 51875,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 76284, 0, 3,
                                                                       75528, 51371, 75612,
                                                                       34343, 34433, 52190,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 76452, 0, 3,
                                                                       75612, 51434, 75696,
                                                                       34433, 34523, 52316,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 76620, 0, 3,
                                                                       75696, 51497, 75780,
                                                                       34523, 34613, 52442,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 76788, 0, 3,
                                                                       75780, 51560, 75864,
                                                                       34613, 34703, 52568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 76956, 0, 3,
                                                                       75864, 51623, 75948,
                                                                       34703, 34793, 52694,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 77124, 0, 3,
                                                                       75948, 51686, 76032,
                                                                       34793, 34883, 52820,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 77292, 0, 3,
                                                                       76032, 51749, 76116,
                                                                       34883, 34973, 52946,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 77460, 0, 3,
                                                                       76116, 51812, 76200,
                                                                       34973, 35063, 53072,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 77628, 0, 3,
                                                                       76284, 52190, 76452,
                                                                       35243, 35393, 53618,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 77908, 0, 3,
                                                                       76452, 52316, 76620,
                                                                       35393, 35543, 53828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 78188, 0, 3,
                                                                       76620, 52442, 76788,
                                                                       35543, 35693, 54038,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 78468, 0, 3,
                                                                       76788, 52568, 76956,
                                                                       35693, 35843, 54248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 78748, 0, 3,
                                                                       76956, 52694, 77124,
                                                                       35843, 35993, 54458,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 79028, 0, 3,
                                                                       77124, 52820, 77292,
                                                                       35993, 36143, 54668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 79308, 0, 3,
                                                                       77292, 52946, 77460,
                                                                       36143, 36293, 54878,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 79588, 0, 3,
                                                                       77628, 53618, 77908,
                                                                       36593, 36818, 55718,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 80008, 0, 3,
                                                                       77908, 53828, 78188,
                                                                       36818, 37043, 56033,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 80428, 0, 3,
                                                                       78188, 54038, 78468,
                                                                       37043, 37268, 56348,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 80848, 0, 3,
                                                                       78468, 54248, 78748,
                                                                       37268, 37493, 56663,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 81268, 0, 3,
                                                                       78748, 54458, 79028,
                                                                       37493, 37718, 56978,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 81688, 0, 3,
                                                                       79028, 54668, 79308,
                                                                       37718, 37943, 57293,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 82108, 0, 3,
                                                                       79588, 55718, 80008,
                                                                       38393, 38708, 58490,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 82696, 0, 3,
                                                                       80008, 56033, 80428,
                                                                       38708, 39023, 58931,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 83284, 0, 3,
                                                                       80428, 56348, 80848,
                                                                       39023, 39338, 59372,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 83872, 0, 3,
                                                                       80848, 56663, 81268,
                                                                       39338, 39653, 59813,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 84460, 0, 3,
                                                                       81268, 56978, 81688,
                                                                       39653, 39968, 60254,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 85048, 0, 3,
                                                                       82108, 58490, 82696,
                                                                       40598, 41018, 61871,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 85832, 0, 3,
                                                                       82696, 58931, 83284,
                                                                       41018, 41438, 62459,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 86616, 0, 3,
                                                                       83284, 59372, 83872,
                                                                       41438, 41858, 63047,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 87400, 0, 3,
                                                                       83872, 59813, 84460,
                                                                       41858, 42278, 63635,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 88184, 0, 3,
                                                                       85048, 61871, 85832,
                                                                       43118, 43658, 65735,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 89192, 0, 3,
                                                                       85832, 62459, 86616,
                                                                       43658, 44198, 66491,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 90200, 0, 3,
                                                                       86616, 63047, 87400,
                                                                       44198, 44738, 67247,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 91208, 0, 3,
                                                                       88184, 65735, 89192,
                                                                       45818, 46493, 69893,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 92468, 0, 3,
                                                                       89192, 66491, 90200,
                                                                       46493, 47168, 70838,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 93728, 0, 3,
                                                                       91208, 69893, 92468,
                                                                       48518, 49343, 74093,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95268, 3, 50993,
                                                                       51014, 75248, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95304, 3, 51014,
                                                                       51035, 75276, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95340, 3, 51035,
                                                                       51056, 75304, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95376, 3, 51056,
                                                                       51077, 75332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95412, 3, 51077,
                                                                       51098, 75360, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95448, 3, 51098,
                                                                       51119, 75388, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95484, 3, 51119,
                                                                       51140, 75416, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95520, 3, 51140,
                                                                       51161, 75444, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95556, 3, 51161,
                                                                       51182, 75472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95592, 3, 51182,
                                                                       51203, 75500, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95628, 0, 3,
                                                                       95268, 75248, 95304,
                                                                       51245, 51308, 75528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95736, 0, 3,
                                                                       95304, 75276, 95340,
                                                                       51308, 51371, 75612,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95844, 0, 3,
                                                                       95340, 75304, 95376,
                                                                       51371, 51434, 75696,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95952, 0, 3,
                                                                       95376, 75332, 95412,
                                                                       51434, 51497, 75780,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 96060, 0, 3,
                                                                       95412, 75360, 95448,
                                                                       51497, 51560, 75864,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 96168, 0, 3,
                                                                       95448, 75388, 95484,
                                                                       51560, 51623, 75948,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 96276, 0, 3,
                                                                       95484, 75416, 95520,
                                                                       51623, 51686, 76032,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 96384, 0, 3,
                                                                       95520, 75444, 95556,
                                                                       51686, 51749, 76116,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 96492, 0, 3,
                                                                       95556, 75472, 95592,
                                                                       51749, 51812, 76200,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 96600, 0, 3,
                                                                       95628, 75528, 95736,
                                                                       51938, 52064, 76284,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 96816, 0, 3,
                                                                       95736, 75612, 95844,
                                                                       52064, 52190, 76452,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97032, 0, 3,
                                                                       95844, 75696, 95952,
                                                                       52190, 52316, 76620,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97248, 0, 3,
                                                                       95952, 75780, 96060,
                                                                       52316, 52442, 76788,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97464, 0, 3,
                                                                       96060, 75864, 96168,
                                                                       52442, 52568, 76956,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97680, 0, 3,
                                                                       96168, 75948, 96276,
                                                                       52568, 52694, 77124,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97896, 0, 3,
                                                                       96276, 76032, 96384,
                                                                       52694, 52820, 77292,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 98112, 0, 3,
                                                                       96384, 76116, 96492,
                                                                       52820, 52946, 77460,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 98328, 0, 3,
                                                                       96600, 76284, 96816,
                                                                       53198, 53408, 77628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 98688, 0, 3,
                                                                       96816, 76452, 97032,
                                                                       53408, 53618, 77908,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 99048, 0, 3,
                                                                       97032, 76620, 97248,
                                                                       53618, 53828, 78188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 99408, 0, 3,
                                                                       97248, 76788, 97464,
                                                                       53828, 54038, 78468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 99768, 0, 3,
                                                                       97464, 76956, 97680,
                                                                       54038, 54248, 78748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 100128, 0, 3,
                                                                       97680, 77124, 97896,
                                                                       54248, 54458, 79028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 100488, 0, 3,
                                                                       97896, 77292, 98112,
                                                                       54458, 54668, 79308,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 100848, 0, 3,
                                                                       98328, 77628, 98688,
                                                                       55088, 55403, 79588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 101388, 0, 3,
                                                                       98688, 77908, 99048,
                                                                       55403, 55718, 80008,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 101928, 0, 3,
                                                                       99048, 78188, 99408,
                                                                       55718, 56033, 80428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 102468, 0, 3,
                                                                       99408, 78468, 99768,
                                                                       56033, 56348, 80848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 103008, 0, 3,
                                                                       99768, 78748, 100128,
                                                                       56348, 56663, 81268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 103548, 0, 3,
                                                                       100128, 79028, 100488,
                                                                       56663, 56978, 81688,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 104088, 0, 3,
                                                                       100848, 79588, 101388,
                                                                       57608, 58049, 82108,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 104844, 0, 3,
                                                                       101388, 80008, 101928,
                                                                       58049, 58490, 82696,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 105600, 0, 3,
                                                                       101928, 80428, 102468,
                                                                       58490, 58931, 83284,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 106356, 0, 3,
                                                                       102468, 80848, 103008,
                                                                       58931, 59372, 83872,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 107112, 0, 3,
                                                                       103008, 81268, 103548,
                                                                       59372, 59813, 84460,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 107868, 0, 3,
                                                                       104088, 82108, 104844,
                                                                       60695, 61283, 85048,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 108876, 0, 3,
                                                                       104844, 82696, 105600,
                                                                       61283, 61871, 85832,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 109884, 0, 3,
                                                                       105600, 83284, 106356,
                                                                       61871, 62459, 86616,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 110892, 0, 3,
                                                                       106356, 83872, 107112,
                                                                       62459, 63047, 87400,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 111900, 0, 3,
                                                                       107868, 85048, 108876,
                                                                       64223, 64979, 88184,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 113196, 0, 3,
                                                                       108876, 85832, 109884,
                                                                       64979, 65735, 89192,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 114492, 0, 3,
                                                                       109884, 86616, 110892,
                                                                       65735, 66491, 90200,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 115788, 0, 3,
                                                                       111900, 88184, 113196,
                                                                       68003, 68948, 91208,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 117408, 0, 3,
                                                                       113196, 89192, 114492,
                                                                       68948, 69893, 92468,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 119028, 0, 3,
                                                                       115788, 91208, 117408,
                                                                       71783, 72938, 93728,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 121008, 104088, 756, ncols);

                    simdfunc::contract_primitives(buffer, 122079, 107868, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 123507, 111900, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 125343, 115788, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 127638, 119028, 1980, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 121764, 121008, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 123087, 122079, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 124803, 123507, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 126963, 125343, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 129618, 127638, 55, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 130443, 121764, 123087, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 131388, 123087, 124803, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 132648, 124803, 126963, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 134268, 126963, 129618, 15,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 136293, 130443, 131388, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 138183, 131388, 132648, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 140703, 132648, 134268, 15,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 143943, 136293, 138183, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 147093, 138183, 140703, 15,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 151293, 143943, 147093, 15,
                                             nmax);

        simdtrf::transform_g_inner(buffer, 156018, 151293, 21, 15, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 156018, 135, nmax);
    }

    for (size_t m = 0; m < 1485; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
