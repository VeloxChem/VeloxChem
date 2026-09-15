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


#include "SimdThreeCenterElectronRepulsionRecIII.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
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
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iii_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iii_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 316030, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2197 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 316030, 194748, 15176, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 18,
                                                             ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 69, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 75, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 78, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 81, 0, 3, 8, 9,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 87, 0, 3, 9, 10,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 93, 0, 3, 10, 11,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 99, 0, 3, 11, 12,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 105, 0, 3, 12, 13,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 111, 0, 3, 13, 14,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 117, 0, 3, 14, 15,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 123, 0, 3, 15, 16,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 129, 0, 3, 16, 17,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 135, 0, 3, 17, 18,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 141, 0, 3, 18, 19,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 147, 0, 3, 19, 20,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 20, 21,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 159, 0, 3, 21, 22,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 165, 0, 3, 22, 23,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 171, 0, 3, 23, 24,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 177, 0, 3, 24, 25,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 27, 30,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 30, 33,
                                                                       87, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 33, 36,
                                                                       93, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 36, 39,
                                                                       99, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 39, 42,
                                                                       105, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 42, 45,
                                                                       111, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 45, 48,
                                                                       117, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 48, 51,
                                                                       123, 129, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 263, 0, 3, 51, 54,
                                                                       129, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 54, 57,
                                                                       135, 141, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 283, 0, 3, 57, 60,
                                                                       141, 147, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 293, 0, 3, 60, 63,
                                                                       147, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 303, 0, 3, 63, 66,
                                                                       153, 159, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 313, 0, 3, 66, 69,
                                                                       159, 165, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 323, 0, 3, 69, 72,
                                                                       165, 171, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 333, 0, 3, 72, 75,
                                                                       171, 177, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 343, 0, 3, 81, 87,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 87, 93,
                                                                       193, 203, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 373, 0, 3, 93, 99,
                                                                       203, 213, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 99,
                                                                       105, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 403, 0, 3, 105,
                                                                       111, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 111,
                                                                       117, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 433, 0, 3, 117,
                                                                       123, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 123,
                                                                       129, 253, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 463, 0, 3, 129,
                                                                       135, 263, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 135,
                                                                       141, 273, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 493, 0, 3, 141,
                                                                       147, 283, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 147,
                                                                       153, 293, 303, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 523, 0, 3, 153,
                                                                       159, 303, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 538, 0, 3, 159,
                                                                       165, 313, 323, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 553, 0, 3, 165,
                                                                       171, 323, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 568, 0, 3, 183,
                                                                       193, 343, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 589, 0, 3, 193,
                                                                       203, 358, 373, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 610, 0, 3, 203,
                                                                       213, 373, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 631, 0, 3, 213,
                                                                       223, 388, 403, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 652, 0, 3, 223,
                                                                       233, 403, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 673, 0, 3, 233,
                                                                       243, 418, 433, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 694, 0, 3, 243,
                                                                       253, 433, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 715, 0, 3, 253,
                                                                       263, 448, 463, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 736, 0, 3, 263,
                                                                       273, 463, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 757, 0, 3, 273,
                                                                       283, 478, 493, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 778, 0, 3, 283,
                                                                       293, 493, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 799, 0, 3, 293,
                                                                       303, 508, 523, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 820, 0, 3, 303,
                                                                       313, 523, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 841, 0, 3, 313,
                                                                       323, 538, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 862, 0, 3, 343,
                                                                       358, 568, 589, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 890, 0, 3, 358,
                                                                       373, 589, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 918, 0, 3, 373,
                                                                       388, 610, 631, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 946, 0, 3, 388,
                                                                       403, 631, 652, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 974, 0, 3, 403,
                                                                       418, 652, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 418,
                                                                       433, 673, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 433,
                                                                       448, 694, 715, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 448,
                                                                       463, 715, 736, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1086, 0, 3, 463,
                                                                       478, 736, 757, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 478,
                                                                       493, 757, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 493,
                                                                       508, 778, 799, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1170, 0, 3, 508,
                                                                       523, 799, 820, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1198, 0, 3, 523,
                                                                       538, 820, 841, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1226, 0, 3, 568,
                                                                       589, 862, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1262, 0, 3, 589,
                                                                       610, 890, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1298, 0, 3, 610,
                                                                       631, 918, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1334, 0, 3, 631,
                                                                       652, 946, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1370, 0, 3, 652,
                                                                       673, 974, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1406, 0, 3, 673,
                                                                       694, 1002, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1442, 0, 3, 694,
                                                                       715, 1030, 1058, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1478, 0, 3, 715,
                                                                       736, 1058, 1086, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1514, 0, 3, 736,
                                                                       757, 1086, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1550, 0, 3, 757,
                                                                       778, 1114, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1586, 0, 3, 778,
                                                                       799, 1142, 1170, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1622, 0, 3, 799,
                                                                       820, 1170, 1198, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1658, 0, 3, 862,
                                                                       890, 1226, 1262, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1703, 0, 3, 890,
                                                                       918, 1262, 1298, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1748, 0, 3, 918,
                                                                       946, 1298, 1334, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1793, 0, 3, 946,
                                                                       974, 1334, 1370, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1838, 0, 3, 974,
                                                                       1002, 1370, 1406, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1883, 0, 3, 1002,
                                                                       1030, 1406, 1442, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1928, 0, 3, 1030,
                                                                       1058, 1442, 1478, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1973, 0, 3, 1058,
                                                                       1086, 1478, 1514, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2018, 0, 3, 1086,
                                                                       1114, 1514, 1550, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2063, 0, 3, 1114,
                                                                       1142, 1550, 1586, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 1142,
                                                                       1170, 1586, 1622, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1226,
                                                                       1262, 1658, 1703, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2208, 0, 3, 1262,
                                                                       1298, 1703, 1748, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2263, 0, 3, 1298,
                                                                       1334, 1748, 1793, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2318, 0, 3, 1334,
                                                                       1370, 1793, 1838, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2373, 0, 3, 1370,
                                                                       1406, 1838, 1883, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2428, 0, 3, 1406,
                                                                       1442, 1883, 1928, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2483, 0, 3, 1442,
                                                                       1478, 1928, 1973, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2538, 0, 3, 1478,
                                                                       1514, 1973, 2018, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2593, 0, 3, 1514,
                                                                       1550, 2018, 2063, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1550,
                                                                       1586, 2063, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2703, 0, 3, 1658,
                                                                       1703, 2153, 2208, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2769, 0, 3, 1703,
                                                                       1748, 2208, 2263, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2835, 0, 3, 1748,
                                                                       1793, 2263, 2318, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2901, 0, 3, 1793,
                                                                       1838, 2318, 2373, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2967, 0, 3, 1838,
                                                                       1883, 2373, 2428, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3033, 0, 3, 1883,
                                                                       1928, 2428, 2483, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3099, 0, 3, 1928,
                                                                       1973, 2483, 2538, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3165, 0, 3, 1973,
                                                                       2018, 2538, 2593, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3231, 0, 3, 2018,
                                                                       2063, 2593, 2648, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3297, 0, 3, 2153,
                                                                       2208, 2703, 2769, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3375, 0, 3, 2208,
                                                                       2263, 2769, 2835, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3453, 0, 3, 2263,
                                                                       2318, 2835, 2901, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3531, 0, 3, 2318,
                                                                       2373, 2901, 2967, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3609, 0, 3, 2373,
                                                                       2428, 2967, 3033, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3687, 0, 3, 2428,
                                                                       2483, 3033, 3099, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3765, 0, 3, 2483,
                                                                       2538, 3099, 3165, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3843, 0, 3, 2538,
                                                                       2593, 3165, 3231, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3921, 0, 3, 2703,
                                                                       2769, 3297, 3375, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4012, 0, 3, 2769,
                                                                       2835, 3375, 3453, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4103, 0, 3, 2835,
                                                                       2901, 3453, 3531, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4194, 0, 3, 2901,
                                                                       2967, 3531, 3609, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4285, 0, 3, 2967,
                                                                       3033, 3609, 3687, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4376, 0, 3, 3033,
                                                                       3099, 3687, 3765, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4467, 0, 3, 3099,
                                                                       3165, 3765, 3843, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4558, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4561, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4564, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4567, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4570, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4573, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4576, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4579, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4582, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4585, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4588, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4591, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4594, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4597, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4600, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4603, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4606, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4609, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4618, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4627, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4636, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4645, 3, 14, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4654, 3, 15, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4663, 3, 16, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4672, 3, 17, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4681, 3, 18, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4690, 3, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4699, 3, 20, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4708, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4717, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4726, 3, 23, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4735, 3, 24, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4744, 3, 25, 78,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4753, 3, 33, 93,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4771, 3, 36, 99,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4789, 3, 39, 105,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4807, 3, 42, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4825, 3, 45, 117,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4843, 3, 48, 123,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4861, 3, 51, 129,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4879, 3, 54, 135,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4897, 3, 57, 141,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4915, 3, 60, 147,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4933, 3, 63, 153,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4951, 3, 66, 159,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4969, 3, 69, 165,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4987, 3, 72, 171,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5005, 3, 75, 177,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5023, 3, 93, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5053, 3, 99, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5083, 3, 105, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5113, 3, 111, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5143, 3, 117, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5173, 3, 123, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5203, 3, 129, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5233, 3, 135, 273,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5263, 3, 141, 283,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5293, 3, 147, 293,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5323, 3, 153, 303,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5353, 3, 159, 313,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5383, 3, 165, 323,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5413, 3, 171, 333,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5443, 3, 203, 373,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5488, 3, 213, 388,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5533, 3, 223, 403,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5578, 3, 233, 418,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5623, 3, 243, 433,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5668, 3, 253, 448,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5713, 3, 263, 463,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5758, 3, 273, 478,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5803, 3, 283, 493,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5848, 3, 293, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5893, 3, 303, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5938, 3, 313, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5983, 3, 323, 553,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6028, 3, 373, 610,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6091, 3, 388, 631,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6154, 3, 403, 652,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6217, 3, 418, 673,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6280, 3, 433, 694,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6343, 3, 448, 715,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6406, 3, 463, 736,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6469, 3, 478, 757,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6532, 3, 493, 778,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6595, 3, 508, 799,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6658, 3, 523, 820,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6721, 3, 538, 841,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6784, 3, 610, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6868, 3, 631, 946,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6952, 3, 652, 974,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7036, 3, 673,
                                                                       1002, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7120, 3, 694,
                                                                       1030, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7204, 3, 715,
                                                                       1058, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7288, 3, 736,
                                                                       1086, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7372, 3, 757,
                                                                       1114, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7456, 3, 778,
                                                                       1142, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7540, 3, 799,
                                                                       1170, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7624, 3, 820,
                                                                       1198, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7708, 3, 918,
                                                                       1298, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7816, 3, 946,
                                                                       1334, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7924, 3, 974,
                                                                       1370, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8032, 3, 1002,
                                                                       1406, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8140, 3, 1030,
                                                                       1442, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8248, 3, 1058,
                                                                       1478, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8356, 3, 1086,
                                                                       1514, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8464, 3, 1114,
                                                                       1550, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8572, 3, 1142,
                                                                       1586, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8680, 3, 1170,
                                                                       1622, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8788, 3, 1298,
                                                                       1748, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8923, 3, 1334,
                                                                       1793, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9058, 3, 1370,
                                                                       1838, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9193, 3, 1406,
                                                                       1883, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9328, 3, 1442,
                                                                       1928, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9463, 3, 1478,
                                                                       1973, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9598, 3, 1514,
                                                                       2018, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9733, 3, 1550,
                                                                       2063, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9868, 3, 1586,
                                                                       2108, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10003, 3, 1748,
                                                                       2263, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10168, 3, 1793,
                                                                       2318, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10333, 3, 1838,
                                                                       2373, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10498, 3, 1883,
                                                                       2428, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10663, 3, 1928,
                                                                       2483, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10828, 3, 1973,
                                                                       2538, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10993, 3, 2018,
                                                                       2593, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11158, 3, 2063,
                                                                       2648, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11323, 3, 2263,
                                                                       2835, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11521, 3, 2318,
                                                                       2901, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11719, 3, 2373,
                                                                       2967, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11917, 3, 2428,
                                                                       3033, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12115, 3, 2483,
                                                                       3099, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12313, 3, 2538,
                                                                       3165, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12511, 3, 2593,
                                                                       3231, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12709, 3, 2835,
                                                                       3453, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12943, 3, 2901,
                                                                       3531, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13177, 3, 2967,
                                                                       3609, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13411, 3, 3033,
                                                                       3687, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13645, 3, 3099,
                                                                       3765, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13879, 3, 3165,
                                                                       3843, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14113, 3, 3453,
                                                                       4103, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14386, 3, 3531,
                                                                       4194, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14659, 3, 3609,
                                                                       4285, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14932, 3, 3687,
                                                                       4376, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 15205, 3, 3765,
                                                                       4467, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15478, 3, 8, 9,
                                                                       4558, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15484, 3, 9, 10,
                                                                       4561, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15490, 3, 10, 11,
                                                                       4564, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15496, 3, 11, 12,
                                                                       4567, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15502, 3, 12, 13,
                                                                       4570, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15508, 3, 13, 14,
                                                                       4573, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15514, 3, 14, 15,
                                                                       4576, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15520, 3, 15, 16,
                                                                       4579, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15526, 3, 16, 17,
                                                                       4582, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15532, 3, 17, 18,
                                                                       4585, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15538, 3, 18, 19,
                                                                       4588, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15544, 3, 19, 20,
                                                                       4591, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15550, 3, 20, 21,
                                                                       4594, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15556, 3, 21, 22,
                                                                       4597, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15562, 3, 22, 23,
                                                                       4600, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15568, 3, 23, 24,
                                                                       4603, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15574, 3, 24, 25,
                                                                       4606, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15580, 0, 3,
                                                                       15478, 4558, 15484, 4609,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15598, 0, 3,
                                                                       15484, 4561, 15490, 4618,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15616, 0, 3,
                                                                       15490, 4564, 15496, 4627,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15634, 0, 3,
                                                                       15496, 4567, 15502, 4636,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15652, 0, 3,
                                                                       15502, 4570, 15508, 4645,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15670, 0, 3,
                                                                       15508, 4573, 15514, 4654,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15688, 0, 3,
                                                                       15514, 4576, 15520, 4663,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15706, 0, 3,
                                                                       15520, 4579, 15526, 4672,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15724, 0, 3,
                                                                       15526, 4582, 15532, 4681,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15742, 0, 3,
                                                                       15532, 4585, 15538, 4690,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15760, 0, 3,
                                                                       15538, 4588, 15544, 4699,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15778, 0, 3,
                                                                       15544, 4591, 15550, 4708,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15796, 0, 3,
                                                                       15550, 4594, 15556, 4717,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15814, 0, 3,
                                                                       15556, 4597, 15562, 4726,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15832, 0, 3,
                                                                       15562, 4600, 15568, 4735,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15850, 0, 3,
                                                                       15568, 4603, 15574, 4744,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15868, 0, 3,
                                                                       15580, 4609, 15598, 81,
                                                                       87, 4753, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15904, 0, 3,
                                                                       15598, 4618, 15616, 87,
                                                                       93, 4771, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15940, 0, 3,
                                                                       15616, 4627, 15634, 93,
                                                                       99, 4789, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15976, 0, 3,
                                                                       15634, 4636, 15652, 99,
                                                                       105, 4807, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16012, 0, 3,
                                                                       15652, 4645, 15670, 105,
                                                                       111, 4825, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16048, 0, 3,
                                                                       15670, 4654, 15688, 111,
                                                                       117, 4843, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16084, 0, 3,
                                                                       15688, 4663, 15706, 117,
                                                                       123, 4861, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16120, 0, 3,
                                                                       15706, 4672, 15724, 123,
                                                                       129, 4879, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16156, 0, 3,
                                                                       15724, 4681, 15742, 129,
                                                                       135, 4897, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16192, 0, 3,
                                                                       15742, 4690, 15760, 135,
                                                                       141, 4915, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16228, 0, 3,
                                                                       15760, 4699, 15778, 141,
                                                                       147, 4933, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16264, 0, 3,
                                                                       15778, 4708, 15796, 147,
                                                                       153, 4951, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16300, 0, 3,
                                                                       15796, 4717, 15814, 153,
                                                                       159, 4969, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16336, 0, 3,
                                                                       15814, 4726, 15832, 159,
                                                                       165, 4987, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16372, 0, 3,
                                                                       15832, 4735, 15850, 165,
                                                                       171, 5005, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16408, 0, 3,
                                                                       15868, 4753, 15904, 183,
                                                                       193, 5023, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16468, 0, 3,
                                                                       15904, 4771, 15940, 193,
                                                                       203, 5053, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16528, 0, 3,
                                                                       15940, 4789, 15976, 203,
                                                                       213, 5083, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16588, 0, 3,
                                                                       15976, 4807, 16012, 213,
                                                                       223, 5113, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16648, 0, 3,
                                                                       16012, 4825, 16048, 223,
                                                                       233, 5143, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16708, 0, 3,
                                                                       16048, 4843, 16084, 233,
                                                                       243, 5173, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16768, 0, 3,
                                                                       16084, 4861, 16120, 243,
                                                                       253, 5203, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16828, 0, 3,
                                                                       16120, 4879, 16156, 253,
                                                                       263, 5233, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16888, 0, 3,
                                                                       16156, 4897, 16192, 263,
                                                                       273, 5263, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16948, 0, 3,
                                                                       16192, 4915, 16228, 273,
                                                                       283, 5293, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17008, 0, 3,
                                                                       16228, 4933, 16264, 283,
                                                                       293, 5323, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17068, 0, 3,
                                                                       16264, 4951, 16300, 293,
                                                                       303, 5353, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17128, 0, 3,
                                                                       16300, 4969, 16336, 303,
                                                                       313, 5383, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17188, 0, 3,
                                                                       16336, 4987, 16372, 313,
                                                                       323, 5413, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17248, 0, 3,
                                                                       16408, 5023, 16468, 343,
                                                                       358, 5443, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17338, 0, 3,
                                                                       16468, 5053, 16528, 358,
                                                                       373, 5488, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17428, 0, 3,
                                                                       16528, 5083, 16588, 373,
                                                                       388, 5533, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17518, 0, 3,
                                                                       16588, 5113, 16648, 388,
                                                                       403, 5578, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17608, 0, 3,
                                                                       16648, 5143, 16708, 403,
                                                                       418, 5623, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17698, 0, 3,
                                                                       16708, 5173, 16768, 418,
                                                                       433, 5668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17788, 0, 3,
                                                                       16768, 5203, 16828, 433,
                                                                       448, 5713, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17878, 0, 3,
                                                                       16828, 5233, 16888, 448,
                                                                       463, 5758, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17968, 0, 3,
                                                                       16888, 5263, 16948, 463,
                                                                       478, 5803, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 18058, 0, 3,
                                                                       16948, 5293, 17008, 478,
                                                                       493, 5848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 18148, 0, 3,
                                                                       17008, 5323, 17068, 493,
                                                                       508, 5893, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 18238, 0, 3,
                                                                       17068, 5353, 17128, 508,
                                                                       523, 5938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 18328, 0, 3,
                                                                       17128, 5383, 17188, 523,
                                                                       538, 5983, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18418, 0, 3,
                                                                       17248, 5443, 17338, 568,
                                                                       589, 6028, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18544, 0, 3,
                                                                       17338, 5488, 17428, 589,
                                                                       610, 6091, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18670, 0, 3,
                                                                       17428, 5533, 17518, 610,
                                                                       631, 6154, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18796, 0, 3,
                                                                       17518, 5578, 17608, 631,
                                                                       652, 6217, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18922, 0, 3,
                                                                       17608, 5623, 17698, 652,
                                                                       673, 6280, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19048, 0, 3,
                                                                       17698, 5668, 17788, 673,
                                                                       694, 6343, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19174, 0, 3,
                                                                       17788, 5713, 17878, 694,
                                                                       715, 6406, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19300, 0, 3,
                                                                       17878, 5758, 17968, 715,
                                                                       736, 6469, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19426, 0, 3,
                                                                       17968, 5803, 18058, 736,
                                                                       757, 6532, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19552, 0, 3,
                                                                       18058, 5848, 18148, 757,
                                                                       778, 6595, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19678, 0, 3,
                                                                       18148, 5893, 18238, 778,
                                                                       799, 6658, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19804, 0, 3,
                                                                       18238, 5938, 18328, 799,
                                                                       820, 6721, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19930, 0, 3,
                                                                       18418, 6028, 18544, 862,
                                                                       890, 6784, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20098, 0, 3,
                                                                       18544, 6091, 18670, 890,
                                                                       918, 6868, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20266, 0, 3,
                                                                       18670, 6154, 18796, 918,
                                                                       946, 6952, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20434, 0, 3,
                                                                       18796, 6217, 18922, 946,
                                                                       974, 7036, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20602, 0, 3,
                                                                       18922, 6280, 19048, 974,
                                                                       1002, 7120, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20770, 0, 3,
                                                                       19048, 6343, 19174, 1002,
                                                                       1030, 7204, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20938, 0, 3,
                                                                       19174, 6406, 19300, 1030,
                                                                       1058, 7288, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21106, 0, 3,
                                                                       19300, 6469, 19426, 1058,
                                                                       1086, 7372, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21274, 0, 3,
                                                                       19426, 6532, 19552, 1086,
                                                                       1114, 7456, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21442, 0, 3,
                                                                       19552, 6595, 19678, 1114,
                                                                       1142, 7540, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21610, 0, 3,
                                                                       19678, 6658, 19804, 1142,
                                                                       1170, 7624, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21778, 0, 3,
                                                                       19930, 6784, 20098, 1226,
                                                                       1262, 7708, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21994, 0, 3,
                                                                       20098, 6868, 20266, 1262,
                                                                       1298, 7816, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22210, 0, 3,
                                                                       20266, 6952, 20434, 1298,
                                                                       1334, 7924, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22426, 0, 3,
                                                                       20434, 7036, 20602, 1334,
                                                                       1370, 8032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22642, 0, 3,
                                                                       20602, 7120, 20770, 1370,
                                                                       1406, 8140, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22858, 0, 3,
                                                                       20770, 7204, 20938, 1406,
                                                                       1442, 8248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23074, 0, 3,
                                                                       20938, 7288, 21106, 1442,
                                                                       1478, 8356, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23290, 0, 3,
                                                                       21106, 7372, 21274, 1478,
                                                                       1514, 8464, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23506, 0, 3,
                                                                       21274, 7456, 21442, 1514,
                                                                       1550, 8572, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23722, 0, 3,
                                                                       21442, 7540, 21610, 1550,
                                                                       1586, 8680, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23938, 0, 3,
                                                                       21778, 7708, 21994, 1658,
                                                                       1703, 8788, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24208, 0, 3,
                                                                       21994, 7816, 22210, 1703,
                                                                       1748, 8923, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24478, 0, 3,
                                                                       22210, 7924, 22426, 1748,
                                                                       1793, 9058, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24748, 0, 3,
                                                                       22426, 8032, 22642, 1793,
                                                                       1838, 9193, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25018, 0, 3,
                                                                       22642, 8140, 22858, 1838,
                                                                       1883, 9328, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25288, 0, 3,
                                                                       22858, 8248, 23074, 1883,
                                                                       1928, 9463, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25558, 0, 3,
                                                                       23074, 8356, 23290, 1928,
                                                                       1973, 9598, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25828, 0, 3,
                                                                       23290, 8464, 23506, 1973,
                                                                       2018, 9733, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26098, 0, 3,
                                                                       23506, 8572, 23722, 2018,
                                                                       2063, 9868, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 26368, 0, 3,
                                                                       23938, 8788, 24208, 2153,
                                                                       2208, 10003, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 26698, 0, 3,
                                                                       24208, 8923, 24478, 2208,
                                                                       2263, 10168, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 27028, 0, 3,
                                                                       24478, 9058, 24748, 2263,
                                                                       2318, 10333, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 27358, 0, 3,
                                                                       24748, 9193, 25018, 2318,
                                                                       2373, 10498, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 27688, 0, 3,
                                                                       25018, 9328, 25288, 2373,
                                                                       2428, 10663, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 28018, 0, 3,
                                                                       25288, 9463, 25558, 2428,
                                                                       2483, 10828, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 28348, 0, 3,
                                                                       25558, 9598, 25828, 2483,
                                                                       2538, 10993, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 28678, 0, 3,
                                                                       25828, 9733, 26098, 2538,
                                                                       2593, 11158, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 29008, 0, 3,
                                                                       26368, 10003, 26698, 2703,
                                                                       2769, 11323, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 29404, 0, 3,
                                                                       26698, 10168, 27028, 2769,
                                                                       2835, 11521, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 29800, 0, 3,
                                                                       27028, 10333, 27358, 2835,
                                                                       2901, 11719, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 30196, 0, 3,
                                                                       27358, 10498, 27688, 2901,
                                                                       2967, 11917, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 30592, 0, 3,
                                                                       27688, 10663, 28018, 2967,
                                                                       3033, 12115, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 30988, 0, 3,
                                                                       28018, 10828, 28348, 3033,
                                                                       3099, 12313, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 31384, 0, 3,
                                                                       28348, 10993, 28678, 3099,
                                                                       3165, 12511, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 31780, 0, 3,
                                                                       29008, 11323, 29404, 3297,
                                                                       3375, 12709, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 32248, 0, 3,
                                                                       29404, 11521, 29800, 3375,
                                                                       3453, 12943, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 32716, 0, 3,
                                                                       29800, 11719, 30196, 3453,
                                                                       3531, 13177, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 33184, 0, 3,
                                                                       30196, 11917, 30592, 3531,
                                                                       3609, 13411, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 33652, 0, 3,
                                                                       30592, 12115, 30988, 3609,
                                                                       3687, 13645, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 34120, 0, 3,
                                                                       30988, 12313, 31384, 3687,
                                                                       3765, 13879, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 34588, 0, 3,
                                                                       31780, 12709, 32248, 3921,
                                                                       4012, 14113, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 35134, 0, 3,
                                                                       32248, 12943, 32716, 4012,
                                                                       4103, 14386, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 35680, 0, 3,
                                                                       32716, 13177, 33184, 4103,
                                                                       4194, 14659, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 36226, 0, 3,
                                                                       33184, 13411, 33652, 4194,
                                                                       4285, 14932, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 36772, 0, 3,
                                                                       33652, 13645, 34120, 4285,
                                                                       4376, 15205, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37318, 3, 4558,
                                                                       4561, 15490, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37328, 3, 4561,
                                                                       4564, 15496, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37338, 3, 4564,
                                                                       4567, 15502, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37348, 3, 4567,
                                                                       4570, 15508, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37358, 3, 4570,
                                                                       4573, 15514, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37368, 3, 4573,
                                                                       4576, 15520, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37378, 3, 4576,
                                                                       4579, 15526, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37388, 3, 4579,
                                                                       4582, 15532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37398, 3, 4582,
                                                                       4585, 15538, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37408, 3, 4585,
                                                                       4588, 15544, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37418, 3, 4588,
                                                                       4591, 15550, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37428, 3, 4591,
                                                                       4594, 15556, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37438, 3, 4594,
                                                                       4597, 15562, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37448, 3, 4597,
                                                                       4600, 15568, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37458, 3, 4600,
                                                                       4603, 15574, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37468, 0, 3,
                                                                       37318, 15490, 37328, 4609,
                                                                       4618, 15616, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37498, 0, 3,
                                                                       37328, 15496, 37338, 4618,
                                                                       4627, 15634, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37528, 0, 3,
                                                                       37338, 15502, 37348, 4627,
                                                                       4636, 15652, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37558, 0, 3,
                                                                       37348, 15508, 37358, 4636,
                                                                       4645, 15670, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37588, 0, 3,
                                                                       37358, 15514, 37368, 4645,
                                                                       4654, 15688, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37618, 0, 3,
                                                                       37368, 15520, 37378, 4654,
                                                                       4663, 15706, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37648, 0, 3,
                                                                       37378, 15526, 37388, 4663,
                                                                       4672, 15724, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37678, 0, 3,
                                                                       37388, 15532, 37398, 4672,
                                                                       4681, 15742, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37708, 0, 3,
                                                                       37398, 15538, 37408, 4681,
                                                                       4690, 15760, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37738, 0, 3,
                                                                       37408, 15544, 37418, 4690,
                                                                       4699, 15778, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37768, 0, 3,
                                                                       37418, 15550, 37428, 4699,
                                                                       4708, 15796, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37798, 0, 3,
                                                                       37428, 15556, 37438, 4708,
                                                                       4717, 15814, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37828, 0, 3,
                                                                       37438, 15562, 37448, 4717,
                                                                       4726, 15832, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37858, 0, 3,
                                                                       37448, 15568, 37458, 4726,
                                                                       4735, 15850, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37888, 0, 3,
                                                                       37468, 15616, 37498, 4753,
                                                                       4771, 15940, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37948, 0, 3,
                                                                       37498, 15634, 37528, 4771,
                                                                       4789, 15976, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38008, 0, 3,
                                                                       37528, 15652, 37558, 4789,
                                                                       4807, 16012, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38068, 0, 3,
                                                                       37558, 15670, 37588, 4807,
                                                                       4825, 16048, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38128, 0, 3,
                                                                       37588, 15688, 37618, 4825,
                                                                       4843, 16084, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38188, 0, 3,
                                                                       37618, 15706, 37648, 4843,
                                                                       4861, 16120, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38248, 0, 3,
                                                                       37648, 15724, 37678, 4861,
                                                                       4879, 16156, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38308, 0, 3,
                                                                       37678, 15742, 37708, 4879,
                                                                       4897, 16192, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38368, 0, 3,
                                                                       37708, 15760, 37738, 4897,
                                                                       4915, 16228, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38428, 0, 3,
                                                                       37738, 15778, 37768, 4915,
                                                                       4933, 16264, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38488, 0, 3,
                                                                       37768, 15796, 37798, 4933,
                                                                       4951, 16300, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38548, 0, 3,
                                                                       37798, 15814, 37828, 4951,
                                                                       4969, 16336, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38608, 0, 3,
                                                                       37828, 15832, 37858, 4969,
                                                                       4987, 16372, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38668, 0, 3,
                                                                       37888, 15940, 37948, 5023,
                                                                       5053, 16528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38768, 0, 3,
                                                                       37948, 15976, 38008, 5053,
                                                                       5083, 16588, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38868, 0, 3,
                                                                       38008, 16012, 38068, 5083,
                                                                       5113, 16648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38968, 0, 3,
                                                                       38068, 16048, 38128, 5113,
                                                                       5143, 16708, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39068, 0, 3,
                                                                       38128, 16084, 38188, 5143,
                                                                       5173, 16768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39168, 0, 3,
                                                                       38188, 16120, 38248, 5173,
                                                                       5203, 16828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39268, 0, 3,
                                                                       38248, 16156, 38308, 5203,
                                                                       5233, 16888, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39368, 0, 3,
                                                                       38308, 16192, 38368, 5233,
                                                                       5263, 16948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39468, 0, 3,
                                                                       38368, 16228, 38428, 5263,
                                                                       5293, 17008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39568, 0, 3,
                                                                       38428, 16264, 38488, 5293,
                                                                       5323, 17068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39668, 0, 3,
                                                                       38488, 16300, 38548, 5323,
                                                                       5353, 17128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39768, 0, 3,
                                                                       38548, 16336, 38608, 5353,
                                                                       5383, 17188, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 39868, 0, 3,
                                                                       38668, 16528, 38768, 5443,
                                                                       5488, 17428, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40018, 0, 3,
                                                                       38768, 16588, 38868, 5488,
                                                                       5533, 17518, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40168, 0, 3,
                                                                       38868, 16648, 38968, 5533,
                                                                       5578, 17608, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40318, 0, 3,
                                                                       38968, 16708, 39068, 5578,
                                                                       5623, 17698, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40468, 0, 3,
                                                                       39068, 16768, 39168, 5623,
                                                                       5668, 17788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40618, 0, 3,
                                                                       39168, 16828, 39268, 5668,
                                                                       5713, 17878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40768, 0, 3,
                                                                       39268, 16888, 39368, 5713,
                                                                       5758, 17968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40918, 0, 3,
                                                                       39368, 16948, 39468, 5758,
                                                                       5803, 18058, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41068, 0, 3,
                                                                       39468, 17008, 39568, 5803,
                                                                       5848, 18148, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41218, 0, 3,
                                                                       39568, 17068, 39668, 5848,
                                                                       5893, 18238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41368, 0, 3,
                                                                       39668, 17128, 39768, 5893,
                                                                       5938, 18328, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 41518, 0, 3,
                                                                       39868, 17428, 40018, 6028,
                                                                       6091, 18670, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 41728, 0, 3,
                                                                       40018, 17518, 40168, 6091,
                                                                       6154, 18796, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 41938, 0, 3,
                                                                       40168, 17608, 40318, 6154,
                                                                       6217, 18922, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 42148, 0, 3,
                                                                       40318, 17698, 40468, 6217,
                                                                       6280, 19048, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 42358, 0, 3,
                                                                       40468, 17788, 40618, 6280,
                                                                       6343, 19174, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 42568, 0, 3,
                                                                       40618, 17878, 40768, 6343,
                                                                       6406, 19300, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 42778, 0, 3,
                                                                       40768, 17968, 40918, 6406,
                                                                       6469, 19426, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 42988, 0, 3,
                                                                       40918, 18058, 41068, 6469,
                                                                       6532, 19552, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43198, 0, 3,
                                                                       41068, 18148, 41218, 6532,
                                                                       6595, 19678, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43408, 0, 3,
                                                                       41218, 18238, 41368, 6595,
                                                                       6658, 19804, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 43618, 0, 3,
                                                                       41518, 18670, 41728, 6784,
                                                                       6868, 20266, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 43898, 0, 3,
                                                                       41728, 18796, 41938, 6868,
                                                                       6952, 20434, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44178, 0, 3,
                                                                       41938, 18922, 42148, 6952,
                                                                       7036, 20602, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44458, 0, 3,
                                                                       42148, 19048, 42358, 7036,
                                                                       7120, 20770, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44738, 0, 3,
                                                                       42358, 19174, 42568, 7120,
                                                                       7204, 20938, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 45018, 0, 3,
                                                                       42568, 19300, 42778, 7204,
                                                                       7288, 21106, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 45298, 0, 3,
                                                                       42778, 19426, 42988, 7288,
                                                                       7372, 21274, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 45578, 0, 3,
                                                                       42988, 19552, 43198, 7372,
                                                                       7456, 21442, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 45858, 0, 3,
                                                                       43198, 19678, 43408, 7456,
                                                                       7540, 21610, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46138, 0, 3,
                                                                       43618, 20266, 43898, 7708,
                                                                       7816, 22210, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46498, 0, 3,
                                                                       43898, 20434, 44178, 7816,
                                                                       7924, 22426, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46858, 0, 3,
                                                                       44178, 20602, 44458, 7924,
                                                                       8032, 22642, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47218, 0, 3,
                                                                       44458, 20770, 44738, 8032,
                                                                       8140, 22858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47578, 0, 3,
                                                                       44738, 20938, 45018, 8140,
                                                                       8248, 23074, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47938, 0, 3,
                                                                       45018, 21106, 45298, 8248,
                                                                       8356, 23290, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 48298, 0, 3,
                                                                       45298, 21274, 45578, 8356,
                                                                       8464, 23506, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 48658, 0, 3,
                                                                       45578, 21442, 45858, 8464,
                                                                       8572, 23722, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 49018, 0, 3,
                                                                       46138, 22210, 46498, 8788,
                                                                       8923, 24478, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 49468, 0, 3,
                                                                       46498, 22426, 46858, 8923,
                                                                       9058, 24748, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 49918, 0, 3,
                                                                       46858, 22642, 47218, 9058,
                                                                       9193, 25018, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 50368, 0, 3,
                                                                       47218, 22858, 47578, 9193,
                                                                       9328, 25288, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 50818, 0, 3,
                                                                       47578, 23074, 47938, 9328,
                                                                       9463, 25558, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 51268, 0, 3,
                                                                       47938, 23290, 48298, 9463,
                                                                       9598, 25828, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 51718, 0, 3,
                                                                       48298, 23506, 48658, 9598,
                                                                       9733, 26098, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 52168, 0, 3,
                                                                       49018, 24478, 49468,
                                                                       10003, 10168, 27028,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 52718, 0, 3,
                                                                       49468, 24748, 49918,
                                                                       10168, 10333, 27358,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 53268, 0, 3,
                                                                       49918, 25018, 50368,
                                                                       10333, 10498, 27688,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 53818, 0, 3,
                                                                       50368, 25288, 50818,
                                                                       10498, 10663, 28018,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 54368, 0, 3,
                                                                       50818, 25558, 51268,
                                                                       10663, 10828, 28348,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 54918, 0, 3,
                                                                       51268, 25828, 51718,
                                                                       10828, 10993, 28678,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 55468, 0, 3,
                                                                       52168, 27028, 52718,
                                                                       11323, 11521, 29800,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 56128, 0, 3,
                                                                       52718, 27358, 53268,
                                                                       11521, 11719, 30196,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 56788, 0, 3,
                                                                       53268, 27688, 53818,
                                                                       11719, 11917, 30592,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 57448, 0, 3,
                                                                       53818, 28018, 54368,
                                                                       11917, 12115, 30988,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 58108, 0, 3,
                                                                       54368, 28348, 54918,
                                                                       12115, 12313, 31384,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 58768, 0, 3,
                                                                       55468, 29800, 56128,
                                                                       12709, 12943, 32716,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 59548, 0, 3,
                                                                       56128, 30196, 56788,
                                                                       12943, 13177, 33184,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 60328, 0, 3,
                                                                       56788, 30592, 57448,
                                                                       13177, 13411, 33652,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 61108, 0, 3,
                                                                       57448, 30988, 58108,
                                                                       13411, 13645, 34120,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 61888, 0, 3,
                                                                       58768, 32716, 59548,
                                                                       14113, 14386, 35680,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 62798, 0, 3,
                                                                       59548, 33184, 60328,
                                                                       14386, 14659, 36226,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 63708, 0, 3,
                                                                       60328, 33652, 61108,
                                                                       14659, 14932, 36772,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64618, 3, 15478,
                                                                       15484, 37318, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64633, 3, 15484,
                                                                       15490, 37328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64648, 3, 15490,
                                                                       15496, 37338, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64663, 3, 15496,
                                                                       15502, 37348, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64678, 3, 15502,
                                                                       15508, 37358, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64693, 3, 15508,
                                                                       15514, 37368, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64708, 3, 15514,
                                                                       15520, 37378, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64723, 3, 15520,
                                                                       15526, 37388, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64738, 3, 15526,
                                                                       15532, 37398, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64753, 3, 15532,
                                                                       15538, 37408, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64768, 3, 15538,
                                                                       15544, 37418, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64783, 3, 15544,
                                                                       15550, 37428, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64798, 3, 15550,
                                                                       15556, 37438, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64813, 3, 15556,
                                                                       15562, 37448, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64828, 3, 15562,
                                                                       15568, 37458, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 64843, 0, 3,
                                                                       64618, 37318, 64633,
                                                                       15580, 15598, 37468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 64888, 0, 3,
                                                                       64633, 37328, 64648,
                                                                       15598, 15616, 37498,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 64933, 0, 3,
                                                                       64648, 37338, 64663,
                                                                       15616, 15634, 37528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 64978, 0, 3,
                                                                       64663, 37348, 64678,
                                                                       15634, 15652, 37558,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65023, 0, 3,
                                                                       64678, 37358, 64693,
                                                                       15652, 15670, 37588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65068, 0, 3,
                                                                       64693, 37368, 64708,
                                                                       15670, 15688, 37618,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65113, 0, 3,
                                                                       64708, 37378, 64723,
                                                                       15688, 15706, 37648,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65158, 0, 3,
                                                                       64723, 37388, 64738,
                                                                       15706, 15724, 37678,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65203, 0, 3,
                                                                       64738, 37398, 64753,
                                                                       15724, 15742, 37708,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65248, 0, 3,
                                                                       64753, 37408, 64768,
                                                                       15742, 15760, 37738,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65293, 0, 3,
                                                                       64768, 37418, 64783,
                                                                       15760, 15778, 37768,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65338, 0, 3,
                                                                       64783, 37428, 64798,
                                                                       15778, 15796, 37798,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65383, 0, 3,
                                                                       64798, 37438, 64813,
                                                                       15796, 15814, 37828,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65428, 0, 3,
                                                                       64813, 37448, 64828,
                                                                       15814, 15832, 37858,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65473, 0, 3,
                                                                       64843, 37468, 64888,
                                                                       15868, 15904, 37888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65563, 0, 3,
                                                                       64888, 37498, 64933,
                                                                       15904, 15940, 37948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65653, 0, 3,
                                                                       64933, 37528, 64978,
                                                                       15940, 15976, 38008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65743, 0, 3,
                                                                       64978, 37558, 65023,
                                                                       15976, 16012, 38068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65833, 0, 3,
                                                                       65023, 37588, 65068,
                                                                       16012, 16048, 38128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65923, 0, 3,
                                                                       65068, 37618, 65113,
                                                                       16048, 16084, 38188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66013, 0, 3,
                                                                       65113, 37648, 65158,
                                                                       16084, 16120, 38248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66103, 0, 3,
                                                                       65158, 37678, 65203,
                                                                       16120, 16156, 38308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66193, 0, 3,
                                                                       65203, 37708, 65248,
                                                                       16156, 16192, 38368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66283, 0, 3,
                                                                       65248, 37738, 65293,
                                                                       16192, 16228, 38428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66373, 0, 3,
                                                                       65293, 37768, 65338,
                                                                       16228, 16264, 38488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66463, 0, 3,
                                                                       65338, 37798, 65383,
                                                                       16264, 16300, 38548,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66553, 0, 3,
                                                                       65383, 37828, 65428,
                                                                       16300, 16336, 38608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 66643, 0, 3,
                                                                       65473, 37888, 65563,
                                                                       16408, 16468, 38668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 66793, 0, 3,
                                                                       65563, 37948, 65653,
                                                                       16468, 16528, 38768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 66943, 0, 3,
                                                                       65653, 38008, 65743,
                                                                       16528, 16588, 38868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67093, 0, 3,
                                                                       65743, 38068, 65833,
                                                                       16588, 16648, 38968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67243, 0, 3,
                                                                       65833, 38128, 65923,
                                                                       16648, 16708, 39068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67393, 0, 3,
                                                                       65923, 38188, 66013,
                                                                       16708, 16768, 39168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67543, 0, 3,
                                                                       66013, 38248, 66103,
                                                                       16768, 16828, 39268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67693, 0, 3,
                                                                       66103, 38308, 66193,
                                                                       16828, 16888, 39368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67843, 0, 3,
                                                                       66193, 38368, 66283,
                                                                       16888, 16948, 39468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67993, 0, 3,
                                                                       66283, 38428, 66373,
                                                                       16948, 17008, 39568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 68143, 0, 3,
                                                                       66373, 38488, 66463,
                                                                       17008, 17068, 39668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 68293, 0, 3,
                                                                       66463, 38548, 66553,
                                                                       17068, 17128, 39768,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 68443, 0, 3,
                                                                       66643, 38668, 66793,
                                                                       17248, 17338, 39868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 68668, 0, 3,
                                                                       66793, 38768, 66943,
                                                                       17338, 17428, 40018,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 68893, 0, 3,
                                                                       66943, 38868, 67093,
                                                                       17428, 17518, 40168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 69118, 0, 3,
                                                                       67093, 38968, 67243,
                                                                       17518, 17608, 40318,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 69343, 0, 3,
                                                                       67243, 39068, 67393,
                                                                       17608, 17698, 40468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 69568, 0, 3,
                                                                       67393, 39168, 67543,
                                                                       17698, 17788, 40618,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 69793, 0, 3,
                                                                       67543, 39268, 67693,
                                                                       17788, 17878, 40768,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 70018, 0, 3,
                                                                       67693, 39368, 67843,
                                                                       17878, 17968, 40918,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 70243, 0, 3,
                                                                       67843, 39468, 67993,
                                                                       17968, 18058, 41068,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 70468, 0, 3,
                                                                       67993, 39568, 68143,
                                                                       18058, 18148, 41218,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 70693, 0, 3,
                                                                       68143, 39668, 68293,
                                                                       18148, 18238, 41368,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 70918, 0, 3,
                                                                       68443, 39868, 68668,
                                                                       18418, 18544, 41518,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 71233, 0, 3,
                                                                       68668, 40018, 68893,
                                                                       18544, 18670, 41728,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 71548, 0, 3,
                                                                       68893, 40168, 69118,
                                                                       18670, 18796, 41938,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 71863, 0, 3,
                                                                       69118, 40318, 69343,
                                                                       18796, 18922, 42148,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 72178, 0, 3,
                                                                       69343, 40468, 69568,
                                                                       18922, 19048, 42358,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 72493, 0, 3,
                                                                       69568, 40618, 69793,
                                                                       19048, 19174, 42568,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 72808, 0, 3,
                                                                       69793, 40768, 70018,
                                                                       19174, 19300, 42778,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 73123, 0, 3,
                                                                       70018, 40918, 70243,
                                                                       19300, 19426, 42988,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 73438, 0, 3,
                                                                       70243, 41068, 70468,
                                                                       19426, 19552, 43198,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 73753, 0, 3,
                                                                       70468, 41218, 70693,
                                                                       19552, 19678, 43408,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 74068, 0, 3,
                                                                       70918, 41518, 71233,
                                                                       19930, 20098, 43618,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 74488, 0, 3,
                                                                       71233, 41728, 71548,
                                                                       20098, 20266, 43898,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 74908, 0, 3,
                                                                       71548, 41938, 71863,
                                                                       20266, 20434, 44178,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 75328, 0, 3,
                                                                       71863, 42148, 72178,
                                                                       20434, 20602, 44458,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 75748, 0, 3,
                                                                       72178, 42358, 72493,
                                                                       20602, 20770, 44738,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 76168, 0, 3,
                                                                       72493, 42568, 72808,
                                                                       20770, 20938, 45018,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 76588, 0, 3,
                                                                       72808, 42778, 73123,
                                                                       20938, 21106, 45298,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 77008, 0, 3,
                                                                       73123, 42988, 73438,
                                                                       21106, 21274, 45578,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 77428, 0, 3,
                                                                       73438, 43198, 73753,
                                                                       21274, 21442, 45858,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 77848, 0, 3,
                                                                       74068, 43618, 74488,
                                                                       21778, 21994, 46138,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 78388, 0, 3,
                                                                       74488, 43898, 74908,
                                                                       21994, 22210, 46498,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 78928, 0, 3,
                                                                       74908, 44178, 75328,
                                                                       22210, 22426, 46858,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 79468, 0, 3,
                                                                       75328, 44458, 75748,
                                                                       22426, 22642, 47218,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 80008, 0, 3,
                                                                       75748, 44738, 76168,
                                                                       22642, 22858, 47578,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 80548, 0, 3,
                                                                       76168, 45018, 76588,
                                                                       22858, 23074, 47938,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 81088, 0, 3,
                                                                       76588, 45298, 77008,
                                                                       23074, 23290, 48298,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 81628, 0, 3,
                                                                       77008, 45578, 77428,
                                                                       23290, 23506, 48658,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 82168, 0, 3,
                                                                       77848, 46138, 78388,
                                                                       23938, 24208, 49018,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 82843, 0, 3,
                                                                       78388, 46498, 78928,
                                                                       24208, 24478, 49468,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 83518, 0, 3,
                                                                       78928, 46858, 79468,
                                                                       24478, 24748, 49918,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 84193, 0, 3,
                                                                       79468, 47218, 80008,
                                                                       24748, 25018, 50368,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 84868, 0, 3,
                                                                       80008, 47578, 80548,
                                                                       25018, 25288, 50818,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 85543, 0, 3,
                                                                       80548, 47938, 81088,
                                                                       25288, 25558, 51268,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 86218, 0, 3,
                                                                       81088, 48298, 81628,
                                                                       25558, 25828, 51718,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 86893, 0, 3,
                                                                       82168, 49018, 82843,
                                                                       26368, 26698, 52168,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 87718, 0, 3,
                                                                       82843, 49468, 83518,
                                                                       26698, 27028, 52718,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 88543, 0, 3,
                                                                       83518, 49918, 84193,
                                                                       27028, 27358, 53268,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 89368, 0, 3,
                                                                       84193, 50368, 84868,
                                                                       27358, 27688, 53818,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 90193, 0, 3,
                                                                       84868, 50818, 85543,
                                                                       27688, 28018, 54368,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 91018, 0, 3,
                                                                       85543, 51268, 86218,
                                                                       28018, 28348, 54918,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 91843, 0, 3,
                                                                       86893, 52168, 87718,
                                                                       29008, 29404, 55468,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 92833, 0, 3,
                                                                       87718, 52718, 88543,
                                                                       29404, 29800, 56128,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 93823, 0, 3,
                                                                       88543, 53268, 89368,
                                                                       29800, 30196, 56788,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 94813, 0, 3,
                                                                       89368, 53818, 90193,
                                                                       30196, 30592, 57448,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 95803, 0, 3,
                                                                       90193, 54368, 91018,
                                                                       30592, 30988, 58108,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 96793, 0, 3,
                                                                       91843, 55468, 92833,
                                                                       31780, 32248, 58768,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 97963, 0, 3,
                                                                       92833, 56128, 93823,
                                                                       32248, 32716, 59548,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 99133, 0, 3,
                                                                       93823, 56788, 94813,
                                                                       32716, 33184, 60328,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 100303, 0, 3,
                                                                       94813, 57448, 95803,
                                                                       33184, 33652, 61108,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 101473, 0, 3,
                                                                       96793, 58768, 97963,
                                                                       34588, 35134, 61888,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 102838, 0, 3,
                                                                       97963, 59548, 99133,
                                                                       35134, 35680, 62798,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 104203, 0, 3,
                                                                       99133, 60328, 100303,
                                                                       35680, 36226, 63708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105568, 3, 37318,
                                                                       37328, 64648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105589, 3, 37328,
                                                                       37338, 64663, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105610, 3, 37338,
                                                                       37348, 64678, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105631, 3, 37348,
                                                                       37358, 64693, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105652, 3, 37358,
                                                                       37368, 64708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105673, 3, 37368,
                                                                       37378, 64723, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105694, 3, 37378,
                                                                       37388, 64738, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105715, 3, 37388,
                                                                       37398, 64753, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105736, 3, 37398,
                                                                       37408, 64768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105757, 3, 37408,
                                                                       37418, 64783, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105778, 3, 37418,
                                                                       37428, 64798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105799, 3, 37428,
                                                                       37438, 64813, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105820, 3, 37438,
                                                                       37448, 64828, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 105841, 0, 3,
                                                                       105568, 64648, 105589,
                                                                       37468, 37498, 64933,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 105904, 0, 3,
                                                                       105589, 64663, 105610,
                                                                       37498, 37528, 64978,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 105967, 0, 3,
                                                                       105610, 64678, 105631,
                                                                       37528, 37558, 65023,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106030, 0, 3,
                                                                       105631, 64693, 105652,
                                                                       37558, 37588, 65068,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106093, 0, 3,
                                                                       105652, 64708, 105673,
                                                                       37588, 37618, 65113,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106156, 0, 3,
                                                                       105673, 64723, 105694,
                                                                       37618, 37648, 65158,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106219, 0, 3,
                                                                       105694, 64738, 105715,
                                                                       37648, 37678, 65203,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106282, 0, 3,
                                                                       105715, 64753, 105736,
                                                                       37678, 37708, 65248,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106345, 0, 3,
                                                                       105736, 64768, 105757,
                                                                       37708, 37738, 65293,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106408, 0, 3,
                                                                       105757, 64783, 105778,
                                                                       37738, 37768, 65338,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106471, 0, 3,
                                                                       105778, 64798, 105799,
                                                                       37768, 37798, 65383,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106534, 0, 3,
                                                                       105799, 64813, 105820,
                                                                       37798, 37828, 65428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 106597, 0, 3,
                                                                       105841, 64933, 105904,
                                                                       37888, 37948, 65653,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 106723, 0, 3,
                                                                       105904, 64978, 105967,
                                                                       37948, 38008, 65743,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 106849, 0, 3,
                                                                       105967, 65023, 106030,
                                                                       38008, 38068, 65833,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 106975, 0, 3,
                                                                       106030, 65068, 106093,
                                                                       38068, 38128, 65923,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107101, 0, 3,
                                                                       106093, 65113, 106156,
                                                                       38128, 38188, 66013,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107227, 0, 3,
                                                                       106156, 65158, 106219,
                                                                       38188, 38248, 66103,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107353, 0, 3,
                                                                       106219, 65203, 106282,
                                                                       38248, 38308, 66193,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107479, 0, 3,
                                                                       106282, 65248, 106345,
                                                                       38308, 38368, 66283,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107605, 0, 3,
                                                                       106345, 65293, 106408,
                                                                       38368, 38428, 66373,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107731, 0, 3,
                                                                       106408, 65338, 106471,
                                                                       38428, 38488, 66463,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107857, 0, 3,
                                                                       106471, 65383, 106534,
                                                                       38488, 38548, 66553,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 107983, 0, 3,
                                                                       106597, 65653, 106723,
                                                                       38668, 38768, 66943,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108193, 0, 3,
                                                                       106723, 65743, 106849,
                                                                       38768, 38868, 67093,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108403, 0, 3,
                                                                       106849, 65833, 106975,
                                                                       38868, 38968, 67243,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108613, 0, 3,
                                                                       106975, 65923, 107101,
                                                                       38968, 39068, 67393,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108823, 0, 3,
                                                                       107101, 66013, 107227,
                                                                       39068, 39168, 67543,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109033, 0, 3,
                                                                       107227, 66103, 107353,
                                                                       39168, 39268, 67693,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109243, 0, 3,
                                                                       107353, 66193, 107479,
                                                                       39268, 39368, 67843,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109453, 0, 3,
                                                                       107479, 66283, 107605,
                                                                       39368, 39468, 67993,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109663, 0, 3,
                                                                       107605, 66373, 107731,
                                                                       39468, 39568, 68143,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109873, 0, 3,
                                                                       107731, 66463, 107857,
                                                                       39568, 39668, 68293,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 110083, 0, 3,
                                                                       107983, 66943, 108193,
                                                                       39868, 40018, 68893,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 110398, 0, 3,
                                                                       108193, 67093, 108403,
                                                                       40018, 40168, 69118,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 110713, 0, 3,
                                                                       108403, 67243, 108613,
                                                                       40168, 40318, 69343,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 111028, 0, 3,
                                                                       108613, 67393, 108823,
                                                                       40318, 40468, 69568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 111343, 0, 3,
                                                                       108823, 67543, 109033,
                                                                       40468, 40618, 69793,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 111658, 0, 3,
                                                                       109033, 67693, 109243,
                                                                       40618, 40768, 70018,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 111973, 0, 3,
                                                                       109243, 67843, 109453,
                                                                       40768, 40918, 70243,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 112288, 0, 3,
                                                                       109453, 67993, 109663,
                                                                       40918, 41068, 70468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 112603, 0, 3,
                                                                       109663, 68143, 109873,
                                                                       41068, 41218, 70693,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 112918, 0, 3,
                                                                       110083, 68893, 110398,
                                                                       41518, 41728, 71548,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 113359, 0, 3,
                                                                       110398, 69118, 110713,
                                                                       41728, 41938, 71863,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 113800, 0, 3,
                                                                       110713, 69343, 111028,
                                                                       41938, 42148, 72178,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 114241, 0, 3,
                                                                       111028, 69568, 111343,
                                                                       42148, 42358, 72493,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 114682, 0, 3,
                                                                       111343, 69793, 111658,
                                                                       42358, 42568, 72808,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 115123, 0, 3,
                                                                       111658, 70018, 111973,
                                                                       42568, 42778, 73123,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 115564, 0, 3,
                                                                       111973, 70243, 112288,
                                                                       42778, 42988, 73438,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 116005, 0, 3,
                                                                       112288, 70468, 112603,
                                                                       42988, 43198, 73753,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 116446, 0, 3,
                                                                       112918, 71548, 113359,
                                                                       43618, 43898, 74908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 117034, 0, 3,
                                                                       113359, 71863, 113800,
                                                                       43898, 44178, 75328,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 117622, 0, 3,
                                                                       113800, 72178, 114241,
                                                                       44178, 44458, 75748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 118210, 0, 3,
                                                                       114241, 72493, 114682,
                                                                       44458, 44738, 76168,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 118798, 0, 3,
                                                                       114682, 72808, 115123,
                                                                       44738, 45018, 76588,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 119386, 0, 3,
                                                                       115123, 73123, 115564,
                                                                       45018, 45298, 77008,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 119974, 0, 3,
                                                                       115564, 73438, 116005,
                                                                       45298, 45578, 77428,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 120562, 0, 3,
                                                                       116446, 74908, 117034,
                                                                       46138, 46498, 78928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 121318, 0, 3,
                                                                       117034, 75328, 117622,
                                                                       46498, 46858, 79468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 122074, 0, 3,
                                                                       117622, 75748, 118210,
                                                                       46858, 47218, 80008,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 122830, 0, 3,
                                                                       118210, 76168, 118798,
                                                                       47218, 47578, 80548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 123586, 0, 3,
                                                                       118798, 76588, 119386,
                                                                       47578, 47938, 81088,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 124342, 0, 3,
                                                                       119386, 77008, 119974,
                                                                       47938, 48298, 81628,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 125098, 0, 3,
                                                                       120562, 78928, 121318,
                                                                       49018, 49468, 83518,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 126043, 0, 3,
                                                                       121318, 79468, 122074,
                                                                       49468, 49918, 84193,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 126988, 0, 3,
                                                                       122074, 80008, 122830,
                                                                       49918, 50368, 84868,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 127933, 0, 3,
                                                                       122830, 80548, 123586,
                                                                       50368, 50818, 85543,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 128878, 0, 3,
                                                                       123586, 81088, 124342,
                                                                       50818, 51268, 86218,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 129823, 0, 3,
                                                                       125098, 83518, 126043,
                                                                       52168, 52718, 88543,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 130978, 0, 3,
                                                                       126043, 84193, 126988,
                                                                       52718, 53268, 89368,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 132133, 0, 3,
                                                                       126988, 84868, 127933,
                                                                       53268, 53818, 90193,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 133288, 0, 3,
                                                                       127933, 85543, 128878,
                                                                       53818, 54368, 91018,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 134443, 0, 3,
                                                                       129823, 88543, 130978,
                                                                       55468, 56128, 93823,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 135829, 0, 3,
                                                                       130978, 89368, 132133,
                                                                       56128, 56788, 94813,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 137215, 0, 3,
                                                                       132133, 90193, 133288,
                                                                       56788, 57448, 95803,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 138601, 0, 3,
                                                                       134443, 93823, 135829,
                                                                       58768, 59548, 99133,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 140239, 0, 3,
                                                                       135829, 94813, 137215,
                                                                       59548, 60328, 100303,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 141877, 0, 3,
                                                                       138601, 99133, 140239,
                                                                       61888, 62798, 104203,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143788, 3, 64618,
                                                                       64633, 105568, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143816, 3, 64633,
                                                                       64648, 105589, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143844, 3, 64648,
                                                                       64663, 105610, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143872, 3, 64663,
                                                                       64678, 105631, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143900, 3, 64678,
                                                                       64693, 105652, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143928, 3, 64693,
                                                                       64708, 105673, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143956, 3, 64708,
                                                                       64723, 105694, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143984, 3, 64723,
                                                                       64738, 105715, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 144012, 3, 64738,
                                                                       64753, 105736, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 144040, 3, 64753,
                                                                       64768, 105757, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 144068, 3, 64768,
                                                                       64783, 105778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 144096, 3, 64783,
                                                                       64798, 105799, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 144124, 3, 64798,
                                                                       64813, 105820, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144152, 0, 3,
                                                                       143788, 105568, 143816,
                                                                       64843, 64888, 105841,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144236, 0, 3,
                                                                       143816, 105589, 143844,
                                                                       64888, 64933, 105904,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144320, 0, 3,
                                                                       143844, 105610, 143872,
                                                                       64933, 64978, 105967,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144404, 0, 3,
                                                                       143872, 105631, 143900,
                                                                       64978, 65023, 106030,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144488, 0, 3,
                                                                       143900, 105652, 143928,
                                                                       65023, 65068, 106093,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144572, 0, 3,
                                                                       143928, 105673, 143956,
                                                                       65068, 65113, 106156,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144656, 0, 3,
                                                                       143956, 105694, 143984,
                                                                       65113, 65158, 106219,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144740, 0, 3,
                                                                       143984, 105715, 144012,
                                                                       65158, 65203, 106282,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144824, 0, 3,
                                                                       144012, 105736, 144040,
                                                                       65203, 65248, 106345,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144908, 0, 3,
                                                                       144040, 105757, 144068,
                                                                       65248, 65293, 106408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144992, 0, 3,
                                                                       144068, 105778, 144096,
                                                                       65293, 65338, 106471,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 145076, 0, 3,
                                                                       144096, 105799, 144124,
                                                                       65338, 65383, 106534,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 145160, 0, 3,
                                                                       144152, 105841, 144236,
                                                                       65473, 65563, 106597,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 145328, 0, 3,
                                                                       144236, 105904, 144320,
                                                                       65563, 65653, 106723,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 145496, 0, 3,
                                                                       144320, 105967, 144404,
                                                                       65653, 65743, 106849,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 145664, 0, 3,
                                                                       144404, 106030, 144488,
                                                                       65743, 65833, 106975,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 145832, 0, 3,
                                                                       144488, 106093, 144572,
                                                                       65833, 65923, 107101,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 146000, 0, 3,
                                                                       144572, 106156, 144656,
                                                                       65923, 66013, 107227,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 146168, 0, 3,
                                                                       144656, 106219, 144740,
                                                                       66013, 66103, 107353,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 146336, 0, 3,
                                                                       144740, 106282, 144824,
                                                                       66103, 66193, 107479,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 146504, 0, 3,
                                                                       144824, 106345, 144908,
                                                                       66193, 66283, 107605,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 146672, 0, 3,
                                                                       144908, 106408, 144992,
                                                                       66283, 66373, 107731,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 146840, 0, 3,
                                                                       144992, 106471, 145076,
                                                                       66373, 66463, 107857,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 147008, 0, 3,
                                                                       145160, 106597, 145328,
                                                                       66643, 66793, 107983,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 147288, 0, 3,
                                                                       145328, 106723, 145496,
                                                                       66793, 66943, 108193,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 147568, 0, 3,
                                                                       145496, 106849, 145664,
                                                                       66943, 67093, 108403,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 147848, 0, 3,
                                                                       145664, 106975, 145832,
                                                                       67093, 67243, 108613,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 148128, 0, 3,
                                                                       145832, 107101, 146000,
                                                                       67243, 67393, 108823,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 148408, 0, 3,
                                                                       146000, 107227, 146168,
                                                                       67393, 67543, 109033,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 148688, 0, 3,
                                                                       146168, 107353, 146336,
                                                                       67543, 67693, 109243,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 148968, 0, 3,
                                                                       146336, 107479, 146504,
                                                                       67693, 67843, 109453,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 149248, 0, 3,
                                                                       146504, 107605, 146672,
                                                                       67843, 67993, 109663,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 149528, 0, 3,
                                                                       146672, 107731, 146840,
                                                                       67993, 68143, 109873,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 149808, 0, 3,
                                                                       147008, 107983, 147288,
                                                                       68443, 68668, 110083,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 150228, 0, 3,
                                                                       147288, 108193, 147568,
                                                                       68668, 68893, 110398,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 150648, 0, 3,
                                                                       147568, 108403, 147848,
                                                                       68893, 69118, 110713,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 151068, 0, 3,
                                                                       147848, 108613, 148128,
                                                                       69118, 69343, 111028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 151488, 0, 3,
                                                                       148128, 108823, 148408,
                                                                       69343, 69568, 111343,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 151908, 0, 3,
                                                                       148408, 109033, 148688,
                                                                       69568, 69793, 111658,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 152328, 0, 3,
                                                                       148688, 109243, 148968,
                                                                       69793, 70018, 111973,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 152748, 0, 3,
                                                                       148968, 109453, 149248,
                                                                       70018, 70243, 112288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 153168, 0, 3,
                                                                       149248, 109663, 149528,
                                                                       70243, 70468, 112603,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 153588, 0, 3,
                                                                       149808, 110083, 150228,
                                                                       70918, 71233, 112918,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 154176, 0, 3,
                                                                       150228, 110398, 150648,
                                                                       71233, 71548, 113359,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 154764, 0, 3,
                                                                       150648, 110713, 151068,
                                                                       71548, 71863, 113800,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 155352, 0, 3,
                                                                       151068, 111028, 151488,
                                                                       71863, 72178, 114241,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 155940, 0, 3,
                                                                       151488, 111343, 151908,
                                                                       72178, 72493, 114682,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 156528, 0, 3,
                                                                       151908, 111658, 152328,
                                                                       72493, 72808, 115123,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 157116, 0, 3,
                                                                       152328, 111973, 152748,
                                                                       72808, 73123, 115564,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 157704, 0, 3,
                                                                       152748, 112288, 153168,
                                                                       73123, 73438, 116005,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 158292, 0, 3,
                                                                       153588, 112918, 154176,
                                                                       74068, 74488, 116446,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 159076, 0, 3,
                                                                       154176, 113359, 154764,
                                                                       74488, 74908, 117034,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 159860, 0, 3,
                                                                       154764, 113800, 155352,
                                                                       74908, 75328, 117622,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 160644, 0, 3,
                                                                       155352, 114241, 155940,
                                                                       75328, 75748, 118210,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 161428, 0, 3,
                                                                       155940, 114682, 156528,
                                                                       75748, 76168, 118798,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 162212, 0, 3,
                                                                       156528, 115123, 157116,
                                                                       76168, 76588, 119386,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 162996, 0, 3,
                                                                       157116, 115564, 157704,
                                                                       76588, 77008, 119974,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 163780, 0, 3,
                                                                       158292, 116446, 159076,
                                                                       77848, 78388, 120562,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 164788, 0, 3,
                                                                       159076, 117034, 159860,
                                                                       78388, 78928, 121318,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 165796, 0, 3,
                                                                       159860, 117622, 160644,
                                                                       78928, 79468, 122074,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 166804, 0, 3,
                                                                       160644, 118210, 161428,
                                                                       79468, 80008, 122830,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 167812, 0, 3,
                                                                       161428, 118798, 162212,
                                                                       80008, 80548, 123586,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 168820, 0, 3,
                                                                       162212, 119386, 162996,
                                                                       80548, 81088, 124342,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 169828, 0, 3,
                                                                       163780, 120562, 164788,
                                                                       82168, 82843, 125098,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 171088, 0, 3,
                                                                       164788, 121318, 165796,
                                                                       82843, 83518, 126043,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 172348, 0, 3,
                                                                       165796, 122074, 166804,
                                                                       83518, 84193, 126988,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 173608, 0, 3,
                                                                       166804, 122830, 167812,
                                                                       84193, 84868, 127933,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 174868, 0, 3,
                                                                       167812, 123586, 168820,
                                                                       84868, 85543, 128878,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 176128, 0, 3,
                                                                       169828, 125098, 171088,
                                                                       86893, 87718, 129823,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 177668, 0, 3,
                                                                       171088, 126043, 172348,
                                                                       87718, 88543, 130978,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 179208, 0, 3,
                                                                       172348, 126988, 173608,
                                                                       88543, 89368, 132133,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 180748, 0, 3,
                                                                       173608, 127933, 174868,
                                                                       89368, 90193, 133288,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 182288, 0, 3,
                                                                       176128, 129823, 177668,
                                                                       91843, 92833, 134443,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 184136, 0, 3,
                                                                       177668, 130978, 179208,
                                                                       92833, 93823, 135829,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 185984, 0, 3,
                                                                       179208, 132133, 180748,
                                                                       93823, 94813, 137215,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 187832, 0, 3,
                                                                       182288, 134443, 184136,
                                                                       96793, 97963, 138601,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 190016, 0, 3,
                                                                       184136, 135829, 185984,
                                                                       97963, 99133, 140239,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 192200, 0, 3,
                                                                       187832, 138601, 190016,
                                                                       101473, 102838, 141877,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 194748, 158292, 784, ncols);

                    simdfunc::contract_primitives(buffer, 195896, 163780, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 197372, 169828, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 199217, 176128, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 201472, 182288, 1848, ncols);

                    simdfunc::contract_primitives(buffer, 204178, 187832, 2184, ncols);

                    simdfunc::contract_primitives(buffer, 207376, 192200, 2548, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 195532, 194748, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 196904, 195896, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 198632, 197372, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 200757, 199217, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 203320, 201472, 66, 1, nmax);

        simdtrf::transform_i_inner(buffer, 206362, 204178, 78, 1, nmax);

        simdtrf::transform_i_inner(buffer, 209924, 207376, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 211107, 195532, 196904, 13,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 212199, 196904, 198632, 13,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 213603, 198632, 200757, 13,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 215358, 200757, 203320, 13,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 217503, 203320, 206362, 13,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 220077, 206362, 209924, 13,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 223119, 211107, 212199, 13,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 225303, 212199, 213603, 13,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 228111, 213603, 215358, 13,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 231621, 215358, 217503, 13,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 235911, 217503, 220077, 13,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 241059, 223119, 225303, 13,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 244699, 225303, 228111, 13,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 249379, 228111, 231621, 13,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 255229, 231621, 235911, 13,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 262379, 241059, 244699, 13,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 267839, 244699, 249379, 13,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 274859, 249379, 255229, 13,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 283634, 262379, 267839, 13,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 291278, 267839, 274859, 13,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 301106, 283634, 291278, 13,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 311298, 301106, 28, 13, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 311298, 169, nmax);
    }

    for (size_t m = 0; m < 2197; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
