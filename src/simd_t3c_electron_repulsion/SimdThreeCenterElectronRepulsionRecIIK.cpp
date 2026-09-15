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


#include "SimdThreeCenterElectronRepulsionRecIIK.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
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
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iik_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iik_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 432622, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2535 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 432622, 291208, 18984, dimensions);

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
                                                        16, 17, 18, 19}, ncols, fj, 6, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4558, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4561, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4564, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4567, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4570, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4573, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4576, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4579, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4582, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4585, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4588, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4591, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4594, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4597, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4600, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4603, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4606, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4609, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4612, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4615, 3, 8, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4624, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4633, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4642, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4651, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4660, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4669, 3, 14, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4678, 3, 15, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4687, 3, 16, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4696, 3, 17, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4705, 3, 18, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4714, 3, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4723, 3, 20, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4732, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4741, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4750, 3, 23, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4759, 3, 24, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4768, 3, 25, 78,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4777, 3, 27, 81,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4795, 3, 30, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4813, 3, 33, 93,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4831, 3, 36, 99,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4849, 3, 39, 105,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4867, 3, 42, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4885, 3, 45, 117,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4903, 3, 48, 123,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4921, 3, 51, 129,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4939, 3, 54, 135,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4957, 3, 57, 141,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4975, 3, 60, 147,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4993, 3, 63, 153,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5011, 3, 66, 159,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5029, 3, 69, 165,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5047, 3, 72, 171,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5065, 3, 75, 177,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5083, 3, 81, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5113, 3, 87, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5143, 3, 93, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5173, 3, 99, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5203, 3, 105, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5233, 3, 111, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5263, 3, 117, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5293, 3, 123, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5323, 3, 129, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5353, 3, 135, 273,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5383, 3, 141, 283,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5413, 3, 147, 293,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5443, 3, 153, 303,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5473, 3, 159, 313,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5503, 3, 165, 323,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5533, 3, 171, 333,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5563, 3, 183, 343,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5608, 3, 193, 358,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5653, 3, 203, 373,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5698, 3, 213, 388,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5743, 3, 223, 403,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5788, 3, 233, 418,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5833, 3, 243, 433,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5878, 3, 253, 448,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5923, 3, 263, 463,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5968, 3, 273, 478,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6013, 3, 283, 493,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6058, 3, 293, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6103, 3, 303, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6148, 3, 313, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6193, 3, 323, 553,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6238, 3, 343, 568,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6301, 3, 358, 589,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6364, 3, 373, 610,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6427, 3, 388, 631,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6490, 3, 403, 652,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6553, 3, 418, 673,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6616, 3, 433, 694,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6679, 3, 448, 715,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6742, 3, 463, 736,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6805, 3, 478, 757,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6868, 3, 493, 778,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6931, 3, 508, 799,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6994, 3, 523, 820,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7057, 3, 538, 841,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7120, 3, 568, 862,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7204, 3, 589, 890,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7288, 3, 610, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7372, 3, 631, 946,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7456, 3, 652, 974,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7540, 3, 673,
                                                                       1002, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7624, 3, 694,
                                                                       1030, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7708, 3, 715,
                                                                       1058, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7792, 3, 736,
                                                                       1086, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7876, 3, 757,
                                                                       1114, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7960, 3, 778,
                                                                       1142, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8044, 3, 799,
                                                                       1170, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8128, 3, 820,
                                                                       1198, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8212, 3, 862,
                                                                       1226, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8320, 3, 890,
                                                                       1262, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8428, 3, 918,
                                                                       1298, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8536, 3, 946,
                                                                       1334, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8644, 3, 974,
                                                                       1370, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8752, 3, 1002,
                                                                       1406, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8860, 3, 1030,
                                                                       1442, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8968, 3, 1058,
                                                                       1478, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9076, 3, 1086,
                                                                       1514, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9184, 3, 1114,
                                                                       1550, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9292, 3, 1142,
                                                                       1586, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9400, 3, 1170,
                                                                       1622, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9508, 3, 1226,
                                                                       1658, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9643, 3, 1262,
                                                                       1703, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9778, 3, 1298,
                                                                       1748, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9913, 3, 1334,
                                                                       1793, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10048, 3, 1370,
                                                                       1838, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10183, 3, 1406,
                                                                       1883, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10318, 3, 1442,
                                                                       1928, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10453, 3, 1478,
                                                                       1973, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10588, 3, 1514,
                                                                       2018, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10723, 3, 1550,
                                                                       2063, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10858, 3, 1586,
                                                                       2108, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10993, 3, 1658,
                                                                       2153, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11158, 3, 1703,
                                                                       2208, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11323, 3, 1748,
                                                                       2263, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11488, 3, 1793,
                                                                       2318, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11653, 3, 1838,
                                                                       2373, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11818, 3, 1883,
                                                                       2428, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11983, 3, 1928,
                                                                       2483, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12148, 3, 1973,
                                                                       2538, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12313, 3, 2018,
                                                                       2593, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12478, 3, 2063,
                                                                       2648, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12643, 3, 2153,
                                                                       2703, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12841, 3, 2208,
                                                                       2769, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13039, 3, 2263,
                                                                       2835, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13237, 3, 2318,
                                                                       2901, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13435, 3, 2373,
                                                                       2967, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13633, 3, 2428,
                                                                       3033, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13831, 3, 2483,
                                                                       3099, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14029, 3, 2538,
                                                                       3165, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14227, 3, 2593,
                                                                       3231, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 14425, 3, 2703,
                                                                       3297, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 14659, 3, 2769,
                                                                       3375, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 14893, 3, 2835,
                                                                       3453, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15127, 3, 2901,
                                                                       3531, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15361, 3, 2967,
                                                                       3609, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15595, 3, 3033,
                                                                       3687, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15829, 3, 3099,
                                                                       3765, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16063, 3, 3165,
                                                                       3843, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 16297, 3, 3297,
                                                                       3921, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 16570, 3, 3375,
                                                                       4012, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 16843, 3, 3453,
                                                                       4103, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17116, 3, 3531,
                                                                       4194, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17389, 3, 3609,
                                                                       4285, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17662, 3, 3687,
                                                                       4376, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17935, 3, 3765,
                                                                       4467, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18208, 3, 8, 9,
                                                                       4564, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18214, 3, 9, 10,
                                                                       4567, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18220, 3, 10, 11,
                                                                       4570, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18226, 3, 11, 12,
                                                                       4573, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18232, 3, 12, 13,
                                                                       4576, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18238, 3, 13, 14,
                                                                       4579, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18244, 3, 14, 15,
                                                                       4582, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18250, 3, 15, 16,
                                                                       4585, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18256, 3, 16, 17,
                                                                       4588, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18262, 3, 17, 18,
                                                                       4591, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18268, 3, 18, 19,
                                                                       4594, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18274, 3, 19, 20,
                                                                       4597, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18280, 3, 20, 21,
                                                                       4600, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18286, 3, 21, 22,
                                                                       4603, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18292, 3, 22, 23,
                                                                       4606, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18298, 3, 23, 24,
                                                                       4609, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18304, 3, 24, 25,
                                                                       4612, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18310, 0, 3,
                                                                       18208, 4564, 18214, 4633,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18328, 0, 3,
                                                                       18214, 4567, 18220, 4642,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18346, 0, 3,
                                                                       18220, 4570, 18226, 4651,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18364, 0, 3,
                                                                       18226, 4573, 18232, 4660,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18382, 0, 3,
                                                                       18232, 4576, 18238, 4669,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18400, 0, 3,
                                                                       18238, 4579, 18244, 4678,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18418, 0, 3,
                                                                       18244, 4582, 18250, 4687,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18436, 0, 3,
                                                                       18250, 4585, 18256, 4696,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18454, 0, 3,
                                                                       18256, 4588, 18262, 4705,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18472, 0, 3,
                                                                       18262, 4591, 18268, 4714,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18490, 0, 3,
                                                                       18268, 4594, 18274, 4723,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18508, 0, 3,
                                                                       18274, 4597, 18280, 4732,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18526, 0, 3,
                                                                       18280, 4600, 18286, 4741,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18544, 0, 3,
                                                                       18286, 4603, 18292, 4750,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18562, 0, 3,
                                                                       18292, 4606, 18298, 4759,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18580, 0, 3,
                                                                       18298, 4609, 18304, 4768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18598, 0, 3,
                                                                       18310, 4633, 18328, 81,
                                                                       87, 4813, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18634, 0, 3,
                                                                       18328, 4642, 18346, 87,
                                                                       93, 4831, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18670, 0, 3,
                                                                       18346, 4651, 18364, 93,
                                                                       99, 4849, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18706, 0, 3,
                                                                       18364, 4660, 18382, 99,
                                                                       105, 4867, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18742, 0, 3,
                                                                       18382, 4669, 18400, 105,
                                                                       111, 4885, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18778, 0, 3,
                                                                       18400, 4678, 18418, 111,
                                                                       117, 4903, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18814, 0, 3,
                                                                       18418, 4687, 18436, 117,
                                                                       123, 4921, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18850, 0, 3,
                                                                       18436, 4696, 18454, 123,
                                                                       129, 4939, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18886, 0, 3,
                                                                       18454, 4705, 18472, 129,
                                                                       135, 4957, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18922, 0, 3,
                                                                       18472, 4714, 18490, 135,
                                                                       141, 4975, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18958, 0, 3,
                                                                       18490, 4723, 18508, 141,
                                                                       147, 4993, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18994, 0, 3,
                                                                       18508, 4732, 18526, 147,
                                                                       153, 5011, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19030, 0, 3,
                                                                       18526, 4741, 18544, 153,
                                                                       159, 5029, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19066, 0, 3,
                                                                       18544, 4750, 18562, 159,
                                                                       165, 5047, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19102, 0, 3,
                                                                       18562, 4759, 18580, 165,
                                                                       171, 5065, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19138, 0, 3,
                                                                       18598, 4813, 18634, 183,
                                                                       193, 5143, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19198, 0, 3,
                                                                       18634, 4831, 18670, 193,
                                                                       203, 5173, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19258, 0, 3,
                                                                       18670, 4849, 18706, 203,
                                                                       213, 5203, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19318, 0, 3,
                                                                       18706, 4867, 18742, 213,
                                                                       223, 5233, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19378, 0, 3,
                                                                       18742, 4885, 18778, 223,
                                                                       233, 5263, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19438, 0, 3,
                                                                       18778, 4903, 18814, 233,
                                                                       243, 5293, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19498, 0, 3,
                                                                       18814, 4921, 18850, 243,
                                                                       253, 5323, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19558, 0, 3,
                                                                       18850, 4939, 18886, 253,
                                                                       263, 5353, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19618, 0, 3,
                                                                       18886, 4957, 18922, 263,
                                                                       273, 5383, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19678, 0, 3,
                                                                       18922, 4975, 18958, 273,
                                                                       283, 5413, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19738, 0, 3,
                                                                       18958, 4993, 18994, 283,
                                                                       293, 5443, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19798, 0, 3,
                                                                       18994, 5011, 19030, 293,
                                                                       303, 5473, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19858, 0, 3,
                                                                       19030, 5029, 19066, 303,
                                                                       313, 5503, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19918, 0, 3,
                                                                       19066, 5047, 19102, 313,
                                                                       323, 5533, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19978, 0, 3,
                                                                       19138, 5143, 19198, 343,
                                                                       358, 5653, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20068, 0, 3,
                                                                       19198, 5173, 19258, 358,
                                                                       373, 5698, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20158, 0, 3,
                                                                       19258, 5203, 19318, 373,
                                                                       388, 5743, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20248, 0, 3,
                                                                       19318, 5233, 19378, 388,
                                                                       403, 5788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20338, 0, 3,
                                                                       19378, 5263, 19438, 403,
                                                                       418, 5833, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20428, 0, 3,
                                                                       19438, 5293, 19498, 418,
                                                                       433, 5878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20518, 0, 3,
                                                                       19498, 5323, 19558, 433,
                                                                       448, 5923, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20608, 0, 3,
                                                                       19558, 5353, 19618, 448,
                                                                       463, 5968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20698, 0, 3,
                                                                       19618, 5383, 19678, 463,
                                                                       478, 6013, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20788, 0, 3,
                                                                       19678, 5413, 19738, 478,
                                                                       493, 6058, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20878, 0, 3,
                                                                       19738, 5443, 19798, 493,
                                                                       508, 6103, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20968, 0, 3,
                                                                       19798, 5473, 19858, 508,
                                                                       523, 6148, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21058, 0, 3,
                                                                       19858, 5503, 19918, 523,
                                                                       538, 6193, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21148, 0, 3,
                                                                       19978, 5653, 20068, 568,
                                                                       589, 6364, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21274, 0, 3,
                                                                       20068, 5698, 20158, 589,
                                                                       610, 6427, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21400, 0, 3,
                                                                       20158, 5743, 20248, 610,
                                                                       631, 6490, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21526, 0, 3,
                                                                       20248, 5788, 20338, 631,
                                                                       652, 6553, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21652, 0, 3,
                                                                       20338, 5833, 20428, 652,
                                                                       673, 6616, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21778, 0, 3,
                                                                       20428, 5878, 20518, 673,
                                                                       694, 6679, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21904, 0, 3,
                                                                       20518, 5923, 20608, 694,
                                                                       715, 6742, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22030, 0, 3,
                                                                       20608, 5968, 20698, 715,
                                                                       736, 6805, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22156, 0, 3,
                                                                       20698, 6013, 20788, 736,
                                                                       757, 6868, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22282, 0, 3,
                                                                       20788, 6058, 20878, 757,
                                                                       778, 6931, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22408, 0, 3,
                                                                       20878, 6103, 20968, 778,
                                                                       799, 6994, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22534, 0, 3,
                                                                       20968, 6148, 21058, 799,
                                                                       820, 7057, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22660, 0, 3,
                                                                       21148, 6364, 21274, 862,
                                                                       890, 7288, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22828, 0, 3,
                                                                       21274, 6427, 21400, 890,
                                                                       918, 7372, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22996, 0, 3,
                                                                       21400, 6490, 21526, 918,
                                                                       946, 7456, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23164, 0, 3,
                                                                       21526, 6553, 21652, 946,
                                                                       974, 7540, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23332, 0, 3,
                                                                       21652, 6616, 21778, 974,
                                                                       1002, 7624, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23500, 0, 3,
                                                                       21778, 6679, 21904, 1002,
                                                                       1030, 7708, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23668, 0, 3,
                                                                       21904, 6742, 22030, 1030,
                                                                       1058, 7792, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23836, 0, 3,
                                                                       22030, 6805, 22156, 1058,
                                                                       1086, 7876, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24004, 0, 3,
                                                                       22156, 6868, 22282, 1086,
                                                                       1114, 7960, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24172, 0, 3,
                                                                       22282, 6931, 22408, 1114,
                                                                       1142, 8044, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24340, 0, 3,
                                                                       22408, 6994, 22534, 1142,
                                                                       1170, 8128, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24508, 0, 3,
                                                                       22660, 7288, 22828, 1226,
                                                                       1262, 8428, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24724, 0, 3,
                                                                       22828, 7372, 22996, 1262,
                                                                       1298, 8536, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24940, 0, 3,
                                                                       22996, 7456, 23164, 1298,
                                                                       1334, 8644, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25156, 0, 3,
                                                                       23164, 7540, 23332, 1334,
                                                                       1370, 8752, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25372, 0, 3,
                                                                       23332, 7624, 23500, 1370,
                                                                       1406, 8860, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25588, 0, 3,
                                                                       23500, 7708, 23668, 1406,
                                                                       1442, 8968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25804, 0, 3,
                                                                       23668, 7792, 23836, 1442,
                                                                       1478, 9076, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26020, 0, 3,
                                                                       23836, 7876, 24004, 1478,
                                                                       1514, 9184, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26236, 0, 3,
                                                                       24004, 7960, 24172, 1514,
                                                                       1550, 9292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26452, 0, 3,
                                                                       24172, 8044, 24340, 1550,
                                                                       1586, 9400, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26668, 0, 3,
                                                                       24508, 8428, 24724, 1658,
                                                                       1703, 9778, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26938, 0, 3,
                                                                       24724, 8536, 24940, 1703,
                                                                       1748, 9913, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27208, 0, 3,
                                                                       24940, 8644, 25156, 1748,
                                                                       1793, 10048, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27478, 0, 3,
                                                                       25156, 8752, 25372, 1793,
                                                                       1838, 10183, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27748, 0, 3,
                                                                       25372, 8860, 25588, 1838,
                                                                       1883, 10318, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28018, 0, 3,
                                                                       25588, 8968, 25804, 1883,
                                                                       1928, 10453, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28288, 0, 3,
                                                                       25804, 9076, 26020, 1928,
                                                                       1973, 10588, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28558, 0, 3,
                                                                       26020, 9184, 26236, 1973,
                                                                       2018, 10723, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28828, 0, 3,
                                                                       26236, 9292, 26452, 2018,
                                                                       2063, 10858, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29098, 0, 3,
                                                                       26668, 9778, 26938, 2153,
                                                                       2208, 11323, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29428, 0, 3,
                                                                       26938, 9913, 27208, 2208,
                                                                       2263, 11488, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29758, 0, 3,
                                                                       27208, 10048, 27478, 2263,
                                                                       2318, 11653, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 30088, 0, 3,
                                                                       27478, 10183, 27748, 2318,
                                                                       2373, 11818, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 30418, 0, 3,
                                                                       27748, 10318, 28018, 2373,
                                                                       2428, 11983, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 30748, 0, 3,
                                                                       28018, 10453, 28288, 2428,
                                                                       2483, 12148, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31078, 0, 3,
                                                                       28288, 10588, 28558, 2483,
                                                                       2538, 12313, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31408, 0, 3,
                                                                       28558, 10723, 28828, 2538,
                                                                       2593, 12478, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 31738, 0, 3,
                                                                       29098, 11323, 29428, 2703,
                                                                       2769, 13039, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 32134, 0, 3,
                                                                       29428, 11488, 29758, 2769,
                                                                       2835, 13237, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 32530, 0, 3,
                                                                       29758, 11653, 30088, 2835,
                                                                       2901, 13435, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 32926, 0, 3,
                                                                       30088, 11818, 30418, 2901,
                                                                       2967, 13633, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 33322, 0, 3,
                                                                       30418, 11983, 30748, 2967,
                                                                       3033, 13831, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 33718, 0, 3,
                                                                       30748, 12148, 31078, 3033,
                                                                       3099, 14029, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 34114, 0, 3,
                                                                       31078, 12313, 31408, 3099,
                                                                       3165, 14227, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 34510, 0, 3,
                                                                       31738, 13039, 32134, 3297,
                                                                       3375, 14893, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 34978, 0, 3,
                                                                       32134, 13237, 32530, 3375,
                                                                       3453, 15127, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 35446, 0, 3,
                                                                       32530, 13435, 32926, 3453,
                                                                       3531, 15361, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 35914, 0, 3,
                                                                       32926, 13633, 33322, 3531,
                                                                       3609, 15595, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 36382, 0, 3,
                                                                       33322, 13831, 33718, 3609,
                                                                       3687, 15829, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 36850, 0, 3,
                                                                       33718, 14029, 34114, 3687,
                                                                       3765, 16063, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 37318, 0, 3,
                                                                       34510, 14893, 34978, 3921,
                                                                       4012, 16843, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 37864, 0, 3,
                                                                       34978, 15127, 35446, 4012,
                                                                       4103, 17116, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 38410, 0, 3,
                                                                       35446, 15361, 35914, 4103,
                                                                       4194, 17389, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 38956, 0, 3,
                                                                       35914, 15595, 36382, 4194,
                                                                       4285, 17662, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 39502, 0, 3,
                                                                       36382, 15829, 36850, 4285,
                                                                       4376, 17935, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40048, 3, 4558,
                                                                       4561, 18208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40058, 3, 4561,
                                                                       4564, 18214, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40068, 3, 4564,
                                                                       4567, 18220, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40078, 3, 4567,
                                                                       4570, 18226, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40088, 3, 4570,
                                                                       4573, 18232, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40098, 3, 4573,
                                                                       4576, 18238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40108, 3, 4576,
                                                                       4579, 18244, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40118, 3, 4579,
                                                                       4582, 18250, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40128, 3, 4582,
                                                                       4585, 18256, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40138, 3, 4585,
                                                                       4588, 18262, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40148, 3, 4588,
                                                                       4591, 18268, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40158, 3, 4591,
                                                                       4594, 18274, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40168, 3, 4594,
                                                                       4597, 18280, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40178, 3, 4597,
                                                                       4600, 18286, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40188, 3, 4600,
                                                                       4603, 18292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40198, 3, 4603,
                                                                       4606, 18298, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40208, 3, 4606,
                                                                       4609, 18304, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40218, 0, 3,
                                                                       40048, 18208, 40058, 4615,
                                                                       4624, 18310, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40248, 0, 3,
                                                                       40058, 18214, 40068, 4624,
                                                                       4633, 18328, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40278, 0, 3,
                                                                       40068, 18220, 40078, 4633,
                                                                       4642, 18346, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40308, 0, 3,
                                                                       40078, 18226, 40088, 4642,
                                                                       4651, 18364, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40338, 0, 3,
                                                                       40088, 18232, 40098, 4651,
                                                                       4660, 18382, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40368, 0, 3,
                                                                       40098, 18238, 40108, 4660,
                                                                       4669, 18400, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40398, 0, 3,
                                                                       40108, 18244, 40118, 4669,
                                                                       4678, 18418, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40428, 0, 3,
                                                                       40118, 18250, 40128, 4678,
                                                                       4687, 18436, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40458, 0, 3,
                                                                       40128, 18256, 40138, 4687,
                                                                       4696, 18454, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40488, 0, 3,
                                                                       40138, 18262, 40148, 4696,
                                                                       4705, 18472, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40518, 0, 3,
                                                                       40148, 18268, 40158, 4705,
                                                                       4714, 18490, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40548, 0, 3,
                                                                       40158, 18274, 40168, 4714,
                                                                       4723, 18508, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40578, 0, 3,
                                                                       40168, 18280, 40178, 4723,
                                                                       4732, 18526, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40608, 0, 3,
                                                                       40178, 18286, 40188, 4732,
                                                                       4741, 18544, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40638, 0, 3,
                                                                       40188, 18292, 40198, 4741,
                                                                       4750, 18562, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40668, 0, 3,
                                                                       40198, 18298, 40208, 4750,
                                                                       4759, 18580, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40698, 0, 3,
                                                                       40218, 18310, 40248, 4777,
                                                                       4795, 18598, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40758, 0, 3,
                                                                       40248, 18328, 40278, 4795,
                                                                       4813, 18634, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40818, 0, 3,
                                                                       40278, 18346, 40308, 4813,
                                                                       4831, 18670, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40878, 0, 3,
                                                                       40308, 18364, 40338, 4831,
                                                                       4849, 18706, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40938, 0, 3,
                                                                       40338, 18382, 40368, 4849,
                                                                       4867, 18742, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40998, 0, 3,
                                                                       40368, 18400, 40398, 4867,
                                                                       4885, 18778, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41058, 0, 3,
                                                                       40398, 18418, 40428, 4885,
                                                                       4903, 18814, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41118, 0, 3,
                                                                       40428, 18436, 40458, 4903,
                                                                       4921, 18850, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41178, 0, 3,
                                                                       40458, 18454, 40488, 4921,
                                                                       4939, 18886, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41238, 0, 3,
                                                                       40488, 18472, 40518, 4939,
                                                                       4957, 18922, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41298, 0, 3,
                                                                       40518, 18490, 40548, 4957,
                                                                       4975, 18958, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41358, 0, 3,
                                                                       40548, 18508, 40578, 4975,
                                                                       4993, 18994, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41418, 0, 3,
                                                                       40578, 18526, 40608, 4993,
                                                                       5011, 19030, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41478, 0, 3,
                                                                       40608, 18544, 40638, 5011,
                                                                       5029, 19066, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41538, 0, 3,
                                                                       40638, 18562, 40668, 5029,
                                                                       5047, 19102, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 41598, 0, 3,
                                                                       40698, 18598, 40758, 5083,
                                                                       5113, 19138, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 41698, 0, 3,
                                                                       40758, 18634, 40818, 5113,
                                                                       5143, 19198, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 41798, 0, 3,
                                                                       40818, 18670, 40878, 5143,
                                                                       5173, 19258, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 41898, 0, 3,
                                                                       40878, 18706, 40938, 5173,
                                                                       5203, 19318, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 41998, 0, 3,
                                                                       40938, 18742, 40998, 5203,
                                                                       5233, 19378, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42098, 0, 3,
                                                                       40998, 18778, 41058, 5233,
                                                                       5263, 19438, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42198, 0, 3,
                                                                       41058, 18814, 41118, 5263,
                                                                       5293, 19498, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42298, 0, 3,
                                                                       41118, 18850, 41178, 5293,
                                                                       5323, 19558, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42398, 0, 3,
                                                                       41178, 18886, 41238, 5323,
                                                                       5353, 19618, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42498, 0, 3,
                                                                       41238, 18922, 41298, 5353,
                                                                       5383, 19678, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42598, 0, 3,
                                                                       41298, 18958, 41358, 5383,
                                                                       5413, 19738, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42698, 0, 3,
                                                                       41358, 18994, 41418, 5413,
                                                                       5443, 19798, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42798, 0, 3,
                                                                       41418, 19030, 41478, 5443,
                                                                       5473, 19858, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42898, 0, 3,
                                                                       41478, 19066, 41538, 5473,
                                                                       5503, 19918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42998, 0, 3,
                                                                       41598, 19138, 41698, 5563,
                                                                       5608, 19978, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43148, 0, 3,
                                                                       41698, 19198, 41798, 5608,
                                                                       5653, 20068, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43298, 0, 3,
                                                                       41798, 19258, 41898, 5653,
                                                                       5698, 20158, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43448, 0, 3,
                                                                       41898, 19318, 41998, 5698,
                                                                       5743, 20248, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43598, 0, 3,
                                                                       41998, 19378, 42098, 5743,
                                                                       5788, 20338, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43748, 0, 3,
                                                                       42098, 19438, 42198, 5788,
                                                                       5833, 20428, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43898, 0, 3,
                                                                       42198, 19498, 42298, 5833,
                                                                       5878, 20518, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44048, 0, 3,
                                                                       42298, 19558, 42398, 5878,
                                                                       5923, 20608, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44198, 0, 3,
                                                                       42398, 19618, 42498, 5923,
                                                                       5968, 20698, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44348, 0, 3,
                                                                       42498, 19678, 42598, 5968,
                                                                       6013, 20788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44498, 0, 3,
                                                                       42598, 19738, 42698, 6013,
                                                                       6058, 20878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44648, 0, 3,
                                                                       42698, 19798, 42798, 6058,
                                                                       6103, 20968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44798, 0, 3,
                                                                       42798, 19858, 42898, 6103,
                                                                       6148, 21058, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 44948, 0, 3,
                                                                       42998, 19978, 43148, 6238,
                                                                       6301, 21148, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45158, 0, 3,
                                                                       43148, 20068, 43298, 6301,
                                                                       6364, 21274, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45368, 0, 3,
                                                                       43298, 20158, 43448, 6364,
                                                                       6427, 21400, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45578, 0, 3,
                                                                       43448, 20248, 43598, 6427,
                                                                       6490, 21526, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45788, 0, 3,
                                                                       43598, 20338, 43748, 6490,
                                                                       6553, 21652, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45998, 0, 3,
                                                                       43748, 20428, 43898, 6553,
                                                                       6616, 21778, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46208, 0, 3,
                                                                       43898, 20518, 44048, 6616,
                                                                       6679, 21904, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46418, 0, 3,
                                                                       44048, 20608, 44198, 6679,
                                                                       6742, 22030, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46628, 0, 3,
                                                                       44198, 20698, 44348, 6742,
                                                                       6805, 22156, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46838, 0, 3,
                                                                       44348, 20788, 44498, 6805,
                                                                       6868, 22282, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 47048, 0, 3,
                                                                       44498, 20878, 44648, 6868,
                                                                       6931, 22408, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 47258, 0, 3,
                                                                       44648, 20968, 44798, 6931,
                                                                       6994, 22534, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 47468, 0, 3,
                                                                       44948, 21148, 45158, 7120,
                                                                       7204, 22660, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 47748, 0, 3,
                                                                       45158, 21274, 45368, 7204,
                                                                       7288, 22828, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48028, 0, 3,
                                                                       45368, 21400, 45578, 7288,
                                                                       7372, 22996, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48308, 0, 3,
                                                                       45578, 21526, 45788, 7372,
                                                                       7456, 23164, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48588, 0, 3,
                                                                       45788, 21652, 45998, 7456,
                                                                       7540, 23332, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48868, 0, 3,
                                                                       45998, 21778, 46208, 7540,
                                                                       7624, 23500, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49148, 0, 3,
                                                                       46208, 21904, 46418, 7624,
                                                                       7708, 23668, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49428, 0, 3,
                                                                       46418, 22030, 46628, 7708,
                                                                       7792, 23836, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49708, 0, 3,
                                                                       46628, 22156, 46838, 7792,
                                                                       7876, 24004, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49988, 0, 3,
                                                                       46838, 22282, 47048, 7876,
                                                                       7960, 24172, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 50268, 0, 3,
                                                                       47048, 22408, 47258, 7960,
                                                                       8044, 24340, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 50548, 0, 3,
                                                                       47468, 22660, 47748, 8212,
                                                                       8320, 24508, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 50908, 0, 3,
                                                                       47748, 22828, 48028, 8320,
                                                                       8428, 24724, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 51268, 0, 3,
                                                                       48028, 22996, 48308, 8428,
                                                                       8536, 24940, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 51628, 0, 3,
                                                                       48308, 23164, 48588, 8536,
                                                                       8644, 25156, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 51988, 0, 3,
                                                                       48588, 23332, 48868, 8644,
                                                                       8752, 25372, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 52348, 0, 3,
                                                                       48868, 23500, 49148, 8752,
                                                                       8860, 25588, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 52708, 0, 3,
                                                                       49148, 23668, 49428, 8860,
                                                                       8968, 25804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53068, 0, 3,
                                                                       49428, 23836, 49708, 8968,
                                                                       9076, 26020, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53428, 0, 3,
                                                                       49708, 24004, 49988, 9076,
                                                                       9184, 26236, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53788, 0, 3,
                                                                       49988, 24172, 50268, 9184,
                                                                       9292, 26452, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 54148, 0, 3,
                                                                       50548, 24508, 50908, 9508,
                                                                       9643, 26668, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 54598, 0, 3,
                                                                       50908, 24724, 51268, 9643,
                                                                       9778, 26938, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 55048, 0, 3,
                                                                       51268, 24940, 51628, 9778,
                                                                       9913, 27208, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 55498, 0, 3,
                                                                       51628, 25156, 51988, 9913,
                                                                       10048, 27478, ncols,
                                                                       gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 55948, 0, 3,
                                                                       51988, 25372, 52348,
                                                                       10048, 10183, 27748,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 56398, 0, 3,
                                                                       52348, 25588, 52708,
                                                                       10183, 10318, 28018,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 56848, 0, 3,
                                                                       52708, 25804, 53068,
                                                                       10318, 10453, 28288,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 57298, 0, 3,
                                                                       53068, 26020, 53428,
                                                                       10453, 10588, 28558,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 57748, 0, 3,
                                                                       53428, 26236, 53788,
                                                                       10588, 10723, 28828,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 58198, 0, 3,
                                                                       54148, 26668, 54598,
                                                                       10993, 11158, 29098,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 58748, 0, 3,
                                                                       54598, 26938, 55048,
                                                                       11158, 11323, 29428,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 59298, 0, 3,
                                                                       55048, 27208, 55498,
                                                                       11323, 11488, 29758,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 59848, 0, 3,
                                                                       55498, 27478, 55948,
                                                                       11488, 11653, 30088,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 60398, 0, 3,
                                                                       55948, 27748, 56398,
                                                                       11653, 11818, 30418,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 60948, 0, 3,
                                                                       56398, 28018, 56848,
                                                                       11818, 11983, 30748,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 61498, 0, 3,
                                                                       56848, 28288, 57298,
                                                                       11983, 12148, 31078,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 62048, 0, 3,
                                                                       57298, 28558, 57748,
                                                                       12148, 12313, 31408,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 62598, 0, 3,
                                                                       58198, 29098, 58748,
                                                                       12643, 12841, 31738,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 63258, 0, 3,
                                                                       58748, 29428, 59298,
                                                                       12841, 13039, 32134,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 63918, 0, 3,
                                                                       59298, 29758, 59848,
                                                                       13039, 13237, 32530,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 64578, 0, 3,
                                                                       59848, 30088, 60398,
                                                                       13237, 13435, 32926,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 65238, 0, 3,
                                                                       60398, 30418, 60948,
                                                                       13435, 13633, 33322,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 65898, 0, 3,
                                                                       60948, 30748, 61498,
                                                                       13633, 13831, 33718,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 66558, 0, 3,
                                                                       61498, 31078, 62048,
                                                                       13831, 14029, 34114,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 67218, 0, 3,
                                                                       62598, 31738, 63258,
                                                                       14425, 14659, 34510,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 67998, 0, 3,
                                                                       63258, 32134, 63918,
                                                                       14659, 14893, 34978,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 68778, 0, 3,
                                                                       63918, 32530, 64578,
                                                                       14893, 15127, 35446,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 69558, 0, 3,
                                                                       64578, 32926, 65238,
                                                                       15127, 15361, 35914,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 70338, 0, 3,
                                                                       65238, 33322, 65898,
                                                                       15361, 15595, 36382,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 71118, 0, 3,
                                                                       65898, 33718, 66558,
                                                                       15595, 15829, 36850,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 71898, 0, 3,
                                                                       67218, 34510, 67998,
                                                                       16297, 16570, 37318,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 72808, 0, 3,
                                                                       67998, 34978, 68778,
                                                                       16570, 16843, 37864,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 73718, 0, 3,
                                                                       68778, 35446, 69558,
                                                                       16843, 17116, 38410,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 74628, 0, 3,
                                                                       69558, 35914, 70338,
                                                                       17116, 17389, 38956,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 75538, 0, 3,
                                                                       70338, 36382, 71118,
                                                                       17389, 17662, 39502,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76448, 3, 18208,
                                                                       18214, 40068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76463, 3, 18214,
                                                                       18220, 40078, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76478, 3, 18220,
                                                                       18226, 40088, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76493, 3, 18226,
                                                                       18232, 40098, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76508, 3, 18232,
                                                                       18238, 40108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76523, 3, 18238,
                                                                       18244, 40118, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76538, 3, 18244,
                                                                       18250, 40128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76553, 3, 18250,
                                                                       18256, 40138, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76568, 3, 18256,
                                                                       18262, 40148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76583, 3, 18262,
                                                                       18268, 40158, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76598, 3, 18268,
                                                                       18274, 40168, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76613, 3, 18274,
                                                                       18280, 40178, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76628, 3, 18280,
                                                                       18286, 40188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76643, 3, 18286,
                                                                       18292, 40198, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76658, 3, 18292,
                                                                       18298, 40208, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76673, 0, 3,
                                                                       76448, 40068, 76463,
                                                                       18310, 18328, 40278,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76718, 0, 3,
                                                                       76463, 40078, 76478,
                                                                       18328, 18346, 40308,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76763, 0, 3,
                                                                       76478, 40088, 76493,
                                                                       18346, 18364, 40338,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76808, 0, 3,
                                                                       76493, 40098, 76508,
                                                                       18364, 18382, 40368,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76853, 0, 3,
                                                                       76508, 40108, 76523,
                                                                       18382, 18400, 40398,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76898, 0, 3,
                                                                       76523, 40118, 76538,
                                                                       18400, 18418, 40428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76943, 0, 3,
                                                                       76538, 40128, 76553,
                                                                       18418, 18436, 40458,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76988, 0, 3,
                                                                       76553, 40138, 76568,
                                                                       18436, 18454, 40488,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77033, 0, 3,
                                                                       76568, 40148, 76583,
                                                                       18454, 18472, 40518,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77078, 0, 3,
                                                                       76583, 40158, 76598,
                                                                       18472, 18490, 40548,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77123, 0, 3,
                                                                       76598, 40168, 76613,
                                                                       18490, 18508, 40578,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77168, 0, 3,
                                                                       76613, 40178, 76628,
                                                                       18508, 18526, 40608,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77213, 0, 3,
                                                                       76628, 40188, 76643,
                                                                       18526, 18544, 40638,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77258, 0, 3,
                                                                       76643, 40198, 76658,
                                                                       18544, 18562, 40668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77303, 0, 3,
                                                                       76673, 40278, 76718,
                                                                       18598, 18634, 40818,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77393, 0, 3,
                                                                       76718, 40308, 76763,
                                                                       18634, 18670, 40878,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77483, 0, 3,
                                                                       76763, 40338, 76808,
                                                                       18670, 18706, 40938,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77573, 0, 3,
                                                                       76808, 40368, 76853,
                                                                       18706, 18742, 40998,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77663, 0, 3,
                                                                       76853, 40398, 76898,
                                                                       18742, 18778, 41058,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77753, 0, 3,
                                                                       76898, 40428, 76943,
                                                                       18778, 18814, 41118,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77843, 0, 3,
                                                                       76943, 40458, 76988,
                                                                       18814, 18850, 41178,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77933, 0, 3,
                                                                       76988, 40488, 77033,
                                                                       18850, 18886, 41238,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78023, 0, 3,
                                                                       77033, 40518, 77078,
                                                                       18886, 18922, 41298,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78113, 0, 3,
                                                                       77078, 40548, 77123,
                                                                       18922, 18958, 41358,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78203, 0, 3,
                                                                       77123, 40578, 77168,
                                                                       18958, 18994, 41418,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78293, 0, 3,
                                                                       77168, 40608, 77213,
                                                                       18994, 19030, 41478,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78383, 0, 3,
                                                                       77213, 40638, 77258,
                                                                       19030, 19066, 41538,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78473, 0, 3,
                                                                       77303, 40818, 77393,
                                                                       19138, 19198, 41798,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78623, 0, 3,
                                                                       77393, 40878, 77483,
                                                                       19198, 19258, 41898,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78773, 0, 3,
                                                                       77483, 40938, 77573,
                                                                       19258, 19318, 41998,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78923, 0, 3,
                                                                       77573, 40998, 77663,
                                                                       19318, 19378, 42098,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79073, 0, 3,
                                                                       77663, 41058, 77753,
                                                                       19378, 19438, 42198,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79223, 0, 3,
                                                                       77753, 41118, 77843,
                                                                       19438, 19498, 42298,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79373, 0, 3,
                                                                       77843, 41178, 77933,
                                                                       19498, 19558, 42398,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79523, 0, 3,
                                                                       77933, 41238, 78023,
                                                                       19558, 19618, 42498,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79673, 0, 3,
                                                                       78023, 41298, 78113,
                                                                       19618, 19678, 42598,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79823, 0, 3,
                                                                       78113, 41358, 78203,
                                                                       19678, 19738, 42698,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79973, 0, 3,
                                                                       78203, 41418, 78293,
                                                                       19738, 19798, 42798,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80123, 0, 3,
                                                                       78293, 41478, 78383,
                                                                       19798, 19858, 42898,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 80273, 0, 3,
                                                                       78473, 41798, 78623,
                                                                       19978, 20068, 43298,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 80498, 0, 3,
                                                                       78623, 41898, 78773,
                                                                       20068, 20158, 43448,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 80723, 0, 3,
                                                                       78773, 41998, 78923,
                                                                       20158, 20248, 43598,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 80948, 0, 3,
                                                                       78923, 42098, 79073,
                                                                       20248, 20338, 43748,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81173, 0, 3,
                                                                       79073, 42198, 79223,
                                                                       20338, 20428, 43898,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81398, 0, 3,
                                                                       79223, 42298, 79373,
                                                                       20428, 20518, 44048,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81623, 0, 3,
                                                                       79373, 42398, 79523,
                                                                       20518, 20608, 44198,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81848, 0, 3,
                                                                       79523, 42498, 79673,
                                                                       20608, 20698, 44348,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82073, 0, 3,
                                                                       79673, 42598, 79823,
                                                                       20698, 20788, 44498,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82298, 0, 3,
                                                                       79823, 42698, 79973,
                                                                       20788, 20878, 44648,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82523, 0, 3,
                                                                       79973, 42798, 80123,
                                                                       20878, 20968, 44798,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 82748, 0, 3,
                                                                       80273, 43298, 80498,
                                                                       21148, 21274, 45368,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 83063, 0, 3,
                                                                       80498, 43448, 80723,
                                                                       21274, 21400, 45578,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 83378, 0, 3,
                                                                       80723, 43598, 80948,
                                                                       21400, 21526, 45788,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 83693, 0, 3,
                                                                       80948, 43748, 81173,
                                                                       21526, 21652, 45998,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 84008, 0, 3,
                                                                       81173, 43898, 81398,
                                                                       21652, 21778, 46208,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 84323, 0, 3,
                                                                       81398, 44048, 81623,
                                                                       21778, 21904, 46418,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 84638, 0, 3,
                                                                       81623, 44198, 81848,
                                                                       21904, 22030, 46628,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 84953, 0, 3,
                                                                       81848, 44348, 82073,
                                                                       22030, 22156, 46838,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 85268, 0, 3,
                                                                       82073, 44498, 82298,
                                                                       22156, 22282, 47048,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 85583, 0, 3,
                                                                       82298, 44648, 82523,
                                                                       22282, 22408, 47258,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 85898, 0, 3,
                                                                       82748, 45368, 83063,
                                                                       22660, 22828, 48028,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 86318, 0, 3,
                                                                       83063, 45578, 83378,
                                                                       22828, 22996, 48308,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 86738, 0, 3,
                                                                       83378, 45788, 83693,
                                                                       22996, 23164, 48588,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 87158, 0, 3,
                                                                       83693, 45998, 84008,
                                                                       23164, 23332, 48868,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 87578, 0, 3,
                                                                       84008, 46208, 84323,
                                                                       23332, 23500, 49148,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 87998, 0, 3,
                                                                       84323, 46418, 84638,
                                                                       23500, 23668, 49428,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 88418, 0, 3,
                                                                       84638, 46628, 84953,
                                                                       23668, 23836, 49708,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 88838, 0, 3,
                                                                       84953, 46838, 85268,
                                                                       23836, 24004, 49988,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 89258, 0, 3,
                                                                       85268, 47048, 85583,
                                                                       24004, 24172, 50268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 89678, 0, 3,
                                                                       85898, 48028, 86318,
                                                                       24508, 24724, 51268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 90218, 0, 3,
                                                                       86318, 48308, 86738,
                                                                       24724, 24940, 51628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 90758, 0, 3,
                                                                       86738, 48588, 87158,
                                                                       24940, 25156, 51988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 91298, 0, 3,
                                                                       87158, 48868, 87578,
                                                                       25156, 25372, 52348,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 91838, 0, 3,
                                                                       87578, 49148, 87998,
                                                                       25372, 25588, 52708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 92378, 0, 3,
                                                                       87998, 49428, 88418,
                                                                       25588, 25804, 53068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 92918, 0, 3,
                                                                       88418, 49708, 88838,
                                                                       25804, 26020, 53428,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 93458, 0, 3,
                                                                       88838, 49988, 89258,
                                                                       26020, 26236, 53788,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 93998, 0, 3,
                                                                       89678, 51268, 90218,
                                                                       26668, 26938, 55048,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 94673, 0, 3,
                                                                       90218, 51628, 90758,
                                                                       26938, 27208, 55498,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 95348, 0, 3,
                                                                       90758, 51988, 91298,
                                                                       27208, 27478, 55948,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 96023, 0, 3,
                                                                       91298, 52348, 91838,
                                                                       27478, 27748, 56398,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 96698, 0, 3,
                                                                       91838, 52708, 92378,
                                                                       27748, 28018, 56848,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 97373, 0, 3,
                                                                       92378, 53068, 92918,
                                                                       28018, 28288, 57298,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 98048, 0, 3,
                                                                       92918, 53428, 93458,
                                                                       28288, 28558, 57748,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 98723, 0, 3,
                                                                       93998, 55048, 94673,
                                                                       29098, 29428, 59298,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 99548, 0, 3,
                                                                       94673, 55498, 95348,
                                                                       29428, 29758, 59848,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 100373, 0, 3,
                                                                       95348, 55948, 96023,
                                                                       29758, 30088, 60398,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 101198, 0, 3,
                                                                       96023, 56398, 96698,
                                                                       30088, 30418, 60948,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 102023, 0, 3,
                                                                       96698, 56848, 97373,
                                                                       30418, 30748, 61498,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 102848, 0, 3,
                                                                       97373, 57298, 98048,
                                                                       30748, 31078, 62048,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 103673, 0, 3,
                                                                       98723, 59298, 99548,
                                                                       31738, 32134, 63918,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 104663, 0, 3,
                                                                       99548, 59848, 100373,
                                                                       32134, 32530, 64578,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 105653, 0, 3,
                                                                       100373, 60398, 101198,
                                                                       32530, 32926, 65238,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 106643, 0, 3,
                                                                       101198, 60948, 102023,
                                                                       32926, 33322, 65898,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 107633, 0, 3,
                                                                       102023, 61498, 102848,
                                                                       33322, 33718, 66558,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 108623, 0, 3,
                                                                       103673, 63918, 104663,
                                                                       34510, 34978, 68778,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 109793, 0, 3,
                                                                       104663, 64578, 105653,
                                                                       34978, 35446, 69558,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 110963, 0, 3,
                                                                       105653, 65238, 106643,
                                                                       35446, 35914, 70338,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 112133, 0, 3,
                                                                       106643, 65898, 107633,
                                                                       35914, 36382, 71118,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 113303, 0, 3,
                                                                       108623, 68778, 109793,
                                                                       37318, 37864, 73718,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 114668, 0, 3,
                                                                       109793, 69558, 110963,
                                                                       37864, 38410, 74628,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 116033, 0, 3,
                                                                       110963, 70338, 112133,
                                                                       38410, 38956, 75538,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117398, 3, 40048,
                                                                       40058, 76448, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117419, 3, 40058,
                                                                       40068, 76463, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117440, 3, 40068,
                                                                       40078, 76478, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117461, 3, 40078,
                                                                       40088, 76493, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117482, 3, 40088,
                                                                       40098, 76508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117503, 3, 40098,
                                                                       40108, 76523, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117524, 3, 40108,
                                                                       40118, 76538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117545, 3, 40118,
                                                                       40128, 76553, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117566, 3, 40128,
                                                                       40138, 76568, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117587, 3, 40138,
                                                                       40148, 76583, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117608, 3, 40148,
                                                                       40158, 76598, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117629, 3, 40158,
                                                                       40168, 76613, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117650, 3, 40168,
                                                                       40178, 76628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117671, 3, 40178,
                                                                       40188, 76643, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117692, 3, 40188,
                                                                       40198, 76658, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 117713, 0, 3,
                                                                       117398, 76448, 117419,
                                                                       40218, 40248, 76673,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 117776, 0, 3,
                                                                       117419, 76463, 117440,
                                                                       40248, 40278, 76718,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 117839, 0, 3,
                                                                       117440, 76478, 117461,
                                                                       40278, 40308, 76763,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 117902, 0, 3,
                                                                       117461, 76493, 117482,
                                                                       40308, 40338, 76808,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 117965, 0, 3,
                                                                       117482, 76508, 117503,
                                                                       40338, 40368, 76853,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118028, 0, 3,
                                                                       117503, 76523, 117524,
                                                                       40368, 40398, 76898,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118091, 0, 3,
                                                                       117524, 76538, 117545,
                                                                       40398, 40428, 76943,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118154, 0, 3,
                                                                       117545, 76553, 117566,
                                                                       40428, 40458, 76988,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118217, 0, 3,
                                                                       117566, 76568, 117587,
                                                                       40458, 40488, 77033,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118280, 0, 3,
                                                                       117587, 76583, 117608,
                                                                       40488, 40518, 77078,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118343, 0, 3,
                                                                       117608, 76598, 117629,
                                                                       40518, 40548, 77123,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118406, 0, 3,
                                                                       117629, 76613, 117650,
                                                                       40548, 40578, 77168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118469, 0, 3,
                                                                       117650, 76628, 117671,
                                                                       40578, 40608, 77213,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118532, 0, 3,
                                                                       117671, 76643, 117692,
                                                                       40608, 40638, 77258,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 118595, 0, 3,
                                                                       117713, 76673, 117776,
                                                                       40698, 40758, 77303,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 118721, 0, 3,
                                                                       117776, 76718, 117839,
                                                                       40758, 40818, 77393,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 118847, 0, 3,
                                                                       117839, 76763, 117902,
                                                                       40818, 40878, 77483,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 118973, 0, 3,
                                                                       117902, 76808, 117965,
                                                                       40878, 40938, 77573,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119099, 0, 3,
                                                                       117965, 76853, 118028,
                                                                       40938, 40998, 77663,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119225, 0, 3,
                                                                       118028, 76898, 118091,
                                                                       40998, 41058, 77753,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119351, 0, 3,
                                                                       118091, 76943, 118154,
                                                                       41058, 41118, 77843,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119477, 0, 3,
                                                                       118154, 76988, 118217,
                                                                       41118, 41178, 77933,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119603, 0, 3,
                                                                       118217, 77033, 118280,
                                                                       41178, 41238, 78023,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119729, 0, 3,
                                                                       118280, 77078, 118343,
                                                                       41238, 41298, 78113,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119855, 0, 3,
                                                                       118343, 77123, 118406,
                                                                       41298, 41358, 78203,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119981, 0, 3,
                                                                       118406, 77168, 118469,
                                                                       41358, 41418, 78293,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 120107, 0, 3,
                                                                       118469, 77213, 118532,
                                                                       41418, 41478, 78383,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 120233, 0, 3,
                                                                       118595, 77303, 118721,
                                                                       41598, 41698, 78473,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 120443, 0, 3,
                                                                       118721, 77393, 118847,
                                                                       41698, 41798, 78623,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 120653, 0, 3,
                                                                       118847, 77483, 118973,
                                                                       41798, 41898, 78773,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 120863, 0, 3,
                                                                       118973, 77573, 119099,
                                                                       41898, 41998, 78923,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 121073, 0, 3,
                                                                       119099, 77663, 119225,
                                                                       41998, 42098, 79073,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 121283, 0, 3,
                                                                       119225, 77753, 119351,
                                                                       42098, 42198, 79223,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 121493, 0, 3,
                                                                       119351, 77843, 119477,
                                                                       42198, 42298, 79373,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 121703, 0, 3,
                                                                       119477, 77933, 119603,
                                                                       42298, 42398, 79523,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 121913, 0, 3,
                                                                       119603, 78023, 119729,
                                                                       42398, 42498, 79673,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 122123, 0, 3,
                                                                       119729, 78113, 119855,
                                                                       42498, 42598, 79823,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 122333, 0, 3,
                                                                       119855, 78203, 119981,
                                                                       42598, 42698, 79973,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 122543, 0, 3,
                                                                       119981, 78293, 120107,
                                                                       42698, 42798, 80123,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 122753, 0, 3,
                                                                       120233, 78473, 120443,
                                                                       42998, 43148, 80273,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 123068, 0, 3,
                                                                       120443, 78623, 120653,
                                                                       43148, 43298, 80498,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 123383, 0, 3,
                                                                       120653, 78773, 120863,
                                                                       43298, 43448, 80723,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 123698, 0, 3,
                                                                       120863, 78923, 121073,
                                                                       43448, 43598, 80948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 124013, 0, 3,
                                                                       121073, 79073, 121283,
                                                                       43598, 43748, 81173,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 124328, 0, 3,
                                                                       121283, 79223, 121493,
                                                                       43748, 43898, 81398,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 124643, 0, 3,
                                                                       121493, 79373, 121703,
                                                                       43898, 44048, 81623,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 124958, 0, 3,
                                                                       121703, 79523, 121913,
                                                                       44048, 44198, 81848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 125273, 0, 3,
                                                                       121913, 79673, 122123,
                                                                       44198, 44348, 82073,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 125588, 0, 3,
                                                                       122123, 79823, 122333,
                                                                       44348, 44498, 82298,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 125903, 0, 3,
                                                                       122333, 79973, 122543,
                                                                       44498, 44648, 82523,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 126218, 0, 3,
                                                                       122753, 80273, 123068,
                                                                       44948, 45158, 82748,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 126659, 0, 3,
                                                                       123068, 80498, 123383,
                                                                       45158, 45368, 83063,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 127100, 0, 3,
                                                                       123383, 80723, 123698,
                                                                       45368, 45578, 83378,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 127541, 0, 3,
                                                                       123698, 80948, 124013,
                                                                       45578, 45788, 83693,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 127982, 0, 3,
                                                                       124013, 81173, 124328,
                                                                       45788, 45998, 84008,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 128423, 0, 3,
                                                                       124328, 81398, 124643,
                                                                       45998, 46208, 84323,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 128864, 0, 3,
                                                                       124643, 81623, 124958,
                                                                       46208, 46418, 84638,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 129305, 0, 3,
                                                                       124958, 81848, 125273,
                                                                       46418, 46628, 84953,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 129746, 0, 3,
                                                                       125273, 82073, 125588,
                                                                       46628, 46838, 85268,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 130187, 0, 3,
                                                                       125588, 82298, 125903,
                                                                       46838, 47048, 85583,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 130628, 0, 3,
                                                                       126218, 82748, 126659,
                                                                       47468, 47748, 85898,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 131216, 0, 3,
                                                                       126659, 83063, 127100,
                                                                       47748, 48028, 86318,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 131804, 0, 3,
                                                                       127100, 83378, 127541,
                                                                       48028, 48308, 86738,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 132392, 0, 3,
                                                                       127541, 83693, 127982,
                                                                       48308, 48588, 87158,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 132980, 0, 3,
                                                                       127982, 84008, 128423,
                                                                       48588, 48868, 87578,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 133568, 0, 3,
                                                                       128423, 84323, 128864,
                                                                       48868, 49148, 87998,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 134156, 0, 3,
                                                                       128864, 84638, 129305,
                                                                       49148, 49428, 88418,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 134744, 0, 3,
                                                                       129305, 84953, 129746,
                                                                       49428, 49708, 88838,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 135332, 0, 3,
                                                                       129746, 85268, 130187,
                                                                       49708, 49988, 89258,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 135920, 0, 3,
                                                                       130628, 85898, 131216,
                                                                       50548, 50908, 89678,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 136676, 0, 3,
                                                                       131216, 86318, 131804,
                                                                       50908, 51268, 90218,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 137432, 0, 3,
                                                                       131804, 86738, 132392,
                                                                       51268, 51628, 90758,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 138188, 0, 3,
                                                                       132392, 87158, 132980,
                                                                       51628, 51988, 91298,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 138944, 0, 3,
                                                                       132980, 87578, 133568,
                                                                       51988, 52348, 91838,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 139700, 0, 3,
                                                                       133568, 87998, 134156,
                                                                       52348, 52708, 92378,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 140456, 0, 3,
                                                                       134156, 88418, 134744,
                                                                       52708, 53068, 92918,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 141212, 0, 3,
                                                                       134744, 88838, 135332,
                                                                       53068, 53428, 93458,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 141968, 0, 3,
                                                                       135920, 89678, 136676,
                                                                       54148, 54598, 93998,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 142913, 0, 3,
                                                                       136676, 90218, 137432,
                                                                       54598, 55048, 94673,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 143858, 0, 3,
                                                                       137432, 90758, 138188,
                                                                       55048, 55498, 95348,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 144803, 0, 3,
                                                                       138188, 91298, 138944,
                                                                       55498, 55948, 96023,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 145748, 0, 3,
                                                                       138944, 91838, 139700,
                                                                       55948, 56398, 96698,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 146693, 0, 3,
                                                                       139700, 92378, 140456,
                                                                       56398, 56848, 97373,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 147638, 0, 3,
                                                                       140456, 92918, 141212,
                                                                       56848, 57298, 98048,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 148583, 0, 3,
                                                                       141968, 93998, 142913,
                                                                       58198, 58748, 98723,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 149738, 0, 3,
                                                                       142913, 94673, 143858,
                                                                       58748, 59298, 99548,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 150893, 0, 3,
                                                                       143858, 95348, 144803,
                                                                       59298, 59848, 100373,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 152048, 0, 3,
                                                                       144803, 96023, 145748,
                                                                       59848, 60398, 101198,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 153203, 0, 3,
                                                                       145748, 96698, 146693,
                                                                       60398, 60948, 102023,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 154358, 0, 3,
                                                                       146693, 97373, 147638,
                                                                       60948, 61498, 102848,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 155513, 0, 3,
                                                                       148583, 98723, 149738,
                                                                       62598, 63258, 103673,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 156899, 0, 3,
                                                                       149738, 99548, 150893,
                                                                       63258, 63918, 104663,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 158285, 0, 3,
                                                                       150893, 100373, 152048,
                                                                       63918, 64578, 105653,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 159671, 0, 3,
                                                                       152048, 101198, 153203,
                                                                       64578, 65238, 106643,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 161057, 0, 3,
                                                                       153203, 102023, 154358,
                                                                       65238, 65898, 107633,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 162443, 0, 3,
                                                                       155513, 103673, 156899,
                                                                       67218, 67998, 108623,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 164081, 0, 3,
                                                                       156899, 104663, 158285,
                                                                       67998, 68778, 109793,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 165719, 0, 3,
                                                                       158285, 105653, 159671,
                                                                       68778, 69558, 110963,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 167357, 0, 3,
                                                                       159671, 106643, 161057,
                                                                       69558, 70338, 112133,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 168995, 0, 3,
                                                                       162443, 108623, 164081,
                                                                       71898, 72808, 113303,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 170906, 0, 3,
                                                                       164081, 109793, 165719,
                                                                       72808, 73718, 114668,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 172817, 0, 3,
                                                                       165719, 110963, 167357,
                                                                       73718, 74628, 116033,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174728, 3, 76448,
                                                                       76463, 117440, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174756, 3, 76463,
                                                                       76478, 117461, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174784, 3, 76478,
                                                                       76493, 117482, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174812, 3, 76493,
                                                                       76508, 117503, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174840, 3, 76508,
                                                                       76523, 117524, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174868, 3, 76523,
                                                                       76538, 117545, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174896, 3, 76538,
                                                                       76553, 117566, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174924, 3, 76553,
                                                                       76568, 117587, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174952, 3, 76568,
                                                                       76583, 117608, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174980, 3, 76583,
                                                                       76598, 117629, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 175008, 3, 76598,
                                                                       76613, 117650, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 175036, 3, 76613,
                                                                       76628, 117671, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 175064, 3, 76628,
                                                                       76643, 117692, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175092, 0, 3,
                                                                       174728, 117440, 174756,
                                                                       76673, 76718, 117839,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175176, 0, 3,
                                                                       174756, 117461, 174784,
                                                                       76718, 76763, 117902,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175260, 0, 3,
                                                                       174784, 117482, 174812,
                                                                       76763, 76808, 117965,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175344, 0, 3,
                                                                       174812, 117503, 174840,
                                                                       76808, 76853, 118028,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175428, 0, 3,
                                                                       174840, 117524, 174868,
                                                                       76853, 76898, 118091,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175512, 0, 3,
                                                                       174868, 117545, 174896,
                                                                       76898, 76943, 118154,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175596, 0, 3,
                                                                       174896, 117566, 174924,
                                                                       76943, 76988, 118217,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175680, 0, 3,
                                                                       174924, 117587, 174952,
                                                                       76988, 77033, 118280,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175764, 0, 3,
                                                                       174952, 117608, 174980,
                                                                       77033, 77078, 118343,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175848, 0, 3,
                                                                       174980, 117629, 175008,
                                                                       77078, 77123, 118406,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175932, 0, 3,
                                                                       175008, 117650, 175036,
                                                                       77123, 77168, 118469,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 176016, 0, 3,
                                                                       175036, 117671, 175064,
                                                                       77168, 77213, 118532,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176100, 0, 3,
                                                                       175092, 117839, 175176,
                                                                       77303, 77393, 118847,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176268, 0, 3,
                                                                       175176, 117902, 175260,
                                                                       77393, 77483, 118973,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176436, 0, 3,
                                                                       175260, 117965, 175344,
                                                                       77483, 77573, 119099,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176604, 0, 3,
                                                                       175344, 118028, 175428,
                                                                       77573, 77663, 119225,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176772, 0, 3,
                                                                       175428, 118091, 175512,
                                                                       77663, 77753, 119351,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176940, 0, 3,
                                                                       175512, 118154, 175596,
                                                                       77753, 77843, 119477,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 177108, 0, 3,
                                                                       175596, 118217, 175680,
                                                                       77843, 77933, 119603,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 177276, 0, 3,
                                                                       175680, 118280, 175764,
                                                                       77933, 78023, 119729,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 177444, 0, 3,
                                                                       175764, 118343, 175848,
                                                                       78023, 78113, 119855,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 177612, 0, 3,
                                                                       175848, 118406, 175932,
                                                                       78113, 78203, 119981,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 177780, 0, 3,
                                                                       175932, 118469, 176016,
                                                                       78203, 78293, 120107,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 177948, 0, 3,
                                                                       176100, 118847, 176268,
                                                                       78473, 78623, 120653,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 178228, 0, 3,
                                                                       176268, 118973, 176436,
                                                                       78623, 78773, 120863,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 178508, 0, 3,
                                                                       176436, 119099, 176604,
                                                                       78773, 78923, 121073,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 178788, 0, 3,
                                                                       176604, 119225, 176772,
                                                                       78923, 79073, 121283,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179068, 0, 3,
                                                                       176772, 119351, 176940,
                                                                       79073, 79223, 121493,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179348, 0, 3,
                                                                       176940, 119477, 177108,
                                                                       79223, 79373, 121703,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179628, 0, 3,
                                                                       177108, 119603, 177276,
                                                                       79373, 79523, 121913,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179908, 0, 3,
                                                                       177276, 119729, 177444,
                                                                       79523, 79673, 122123,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 180188, 0, 3,
                                                                       177444, 119855, 177612,
                                                                       79673, 79823, 122333,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 180468, 0, 3,
                                                                       177612, 119981, 177780,
                                                                       79823, 79973, 122543,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 180748, 0, 3,
                                                                       177948, 120653, 178228,
                                                                       80273, 80498, 123383,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 181168, 0, 3,
                                                                       178228, 120863, 178508,
                                                                       80498, 80723, 123698,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 181588, 0, 3,
                                                                       178508, 121073, 178788,
                                                                       80723, 80948, 124013,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 182008, 0, 3,
                                                                       178788, 121283, 179068,
                                                                       80948, 81173, 124328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 182428, 0, 3,
                                                                       179068, 121493, 179348,
                                                                       81173, 81398, 124643,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 182848, 0, 3,
                                                                       179348, 121703, 179628,
                                                                       81398, 81623, 124958,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 183268, 0, 3,
                                                                       179628, 121913, 179908,
                                                                       81623, 81848, 125273,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 183688, 0, 3,
                                                                       179908, 122123, 180188,
                                                                       81848, 82073, 125588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 184108, 0, 3,
                                                                       180188, 122333, 180468,
                                                                       82073, 82298, 125903,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 184528, 0, 3,
                                                                       180748, 123383, 181168,
                                                                       82748, 83063, 127100,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 185116, 0, 3,
                                                                       181168, 123698, 181588,
                                                                       83063, 83378, 127541,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 185704, 0, 3,
                                                                       181588, 124013, 182008,
                                                                       83378, 83693, 127982,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 186292, 0, 3,
                                                                       182008, 124328, 182428,
                                                                       83693, 84008, 128423,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 186880, 0, 3,
                                                                       182428, 124643, 182848,
                                                                       84008, 84323, 128864,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 187468, 0, 3,
                                                                       182848, 124958, 183268,
                                                                       84323, 84638, 129305,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 188056, 0, 3,
                                                                       183268, 125273, 183688,
                                                                       84638, 84953, 129746,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 188644, 0, 3,
                                                                       183688, 125588, 184108,
                                                                       84953, 85268, 130187,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 189232, 0, 3,
                                                                       184528, 127100, 185116,
                                                                       85898, 86318, 131804,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 190016, 0, 3,
                                                                       185116, 127541, 185704,
                                                                       86318, 86738, 132392,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 190800, 0, 3,
                                                                       185704, 127982, 186292,
                                                                       86738, 87158, 132980,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 191584, 0, 3,
                                                                       186292, 128423, 186880,
                                                                       87158, 87578, 133568,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 192368, 0, 3,
                                                                       186880, 128864, 187468,
                                                                       87578, 87998, 134156,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 193152, 0, 3,
                                                                       187468, 129305, 188056,
                                                                       87998, 88418, 134744,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 193936, 0, 3,
                                                                       188056, 129746, 188644,
                                                                       88418, 88838, 135332,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 194720, 0, 3,
                                                                       189232, 131804, 190016,
                                                                       89678, 90218, 137432,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 195728, 0, 3,
                                                                       190016, 132392, 190800,
                                                                       90218, 90758, 138188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 196736, 0, 3,
                                                                       190800, 132980, 191584,
                                                                       90758, 91298, 138944,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 197744, 0, 3,
                                                                       191584, 133568, 192368,
                                                                       91298, 91838, 139700,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 198752, 0, 3,
                                                                       192368, 134156, 193152,
                                                                       91838, 92378, 140456,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 199760, 0, 3,
                                                                       193152, 134744, 193936,
                                                                       92378, 92918, 141212,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 200768, 0, 3,
                                                                       194720, 137432, 195728,
                                                                       93998, 94673, 143858,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 202028, 0, 3,
                                                                       195728, 138188, 196736,
                                                                       94673, 95348, 144803,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 203288, 0, 3,
                                                                       196736, 138944, 197744,
                                                                       95348, 96023, 145748,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 204548, 0, 3,
                                                                       197744, 139700, 198752,
                                                                       96023, 96698, 146693,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 205808, 0, 3,
                                                                       198752, 140456, 199760,
                                                                       96698, 97373, 147638,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 207068, 0, 3,
                                                                       200768, 143858, 202028,
                                                                       98723, 99548, 150893,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 208608, 0, 3,
                                                                       202028, 144803, 203288,
                                                                       99548, 100373, 152048,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 210148, 0, 3,
                                                                       203288, 145748, 204548,
                                                                       100373, 101198, 153203,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 211688, 0, 3,
                                                                       204548, 146693, 205808,
                                                                       101198, 102023, 154358,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 213228, 0, 3,
                                                                       207068, 150893, 208608,
                                                                       103673, 104663, 158285,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 215076, 0, 3,
                                                                       208608, 152048, 210148,
                                                                       104663, 105653, 159671,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 216924, 0, 3,
                                                                       210148, 153203, 211688,
                                                                       105653, 106643, 161057,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 218772, 0, 3,
                                                                       213228, 158285, 215076,
                                                                       108623, 109793, 165719,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 220956, 0, 3,
                                                                       215076, 159671, 216924,
                                                                       109793, 110963, 167357,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 223140, 0, 3,
                                                                       218772, 165719, 220956,
                                                                       113303, 114668, 172817,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225688, 3, 117398,
                                                                       117419, 174728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225724, 3, 117419,
                                                                       117440, 174756, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225760, 3, 117440,
                                                                       117461, 174784, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225796, 3, 117461,
                                                                       117482, 174812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225832, 3, 117482,
                                                                       117503, 174840, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225868, 3, 117503,
                                                                       117524, 174868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225904, 3, 117524,
                                                                       117545, 174896, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225940, 3, 117545,
                                                                       117566, 174924, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225976, 3, 117566,
                                                                       117587, 174952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 226012, 3, 117587,
                                                                       117608, 174980, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 226048, 3, 117608,
                                                                       117629, 175008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 226084, 3, 117629,
                                                                       117650, 175036, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 226120, 3, 117650,
                                                                       117671, 175064, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226156, 0, 3,
                                                                       225688, 174728, 225724,
                                                                       117713, 117776, 175092,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226264, 0, 3,
                                                                       225724, 174756, 225760,
                                                                       117776, 117839, 175176,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226372, 0, 3,
                                                                       225760, 174784, 225796,
                                                                       117839, 117902, 175260,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226480, 0, 3,
                                                                       225796, 174812, 225832,
                                                                       117902, 117965, 175344,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226588, 0, 3,
                                                                       225832, 174840, 225868,
                                                                       117965, 118028, 175428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226696, 0, 3,
                                                                       225868, 174868, 225904,
                                                                       118028, 118091, 175512,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226804, 0, 3,
                                                                       225904, 174896, 225940,
                                                                       118091, 118154, 175596,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226912, 0, 3,
                                                                       225940, 174924, 225976,
                                                                       118154, 118217, 175680,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 227020, 0, 3,
                                                                       225976, 174952, 226012,
                                                                       118217, 118280, 175764,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 227128, 0, 3,
                                                                       226012, 174980, 226048,
                                                                       118280, 118343, 175848,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 227236, 0, 3,
                                                                       226048, 175008, 226084,
                                                                       118343, 118406, 175932,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 227344, 0, 3,
                                                                       226084, 175036, 226120,
                                                                       118406, 118469, 176016,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 227452, 0, 3,
                                                                       226156, 175092, 226264,
                                                                       118595, 118721, 176100,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 227668, 0, 3,
                                                                       226264, 175176, 226372,
                                                                       118721, 118847, 176268,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 227884, 0, 3,
                                                                       226372, 175260, 226480,
                                                                       118847, 118973, 176436,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 228100, 0, 3,
                                                                       226480, 175344, 226588,
                                                                       118973, 119099, 176604,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 228316, 0, 3,
                                                                       226588, 175428, 226696,
                                                                       119099, 119225, 176772,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 228532, 0, 3,
                                                                       226696, 175512, 226804,
                                                                       119225, 119351, 176940,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 228748, 0, 3,
                                                                       226804, 175596, 226912,
                                                                       119351, 119477, 177108,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 228964, 0, 3,
                                                                       226912, 175680, 227020,
                                                                       119477, 119603, 177276,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 229180, 0, 3,
                                                                       227020, 175764, 227128,
                                                                       119603, 119729, 177444,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 229396, 0, 3,
                                                                       227128, 175848, 227236,
                                                                       119729, 119855, 177612,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 229612, 0, 3,
                                                                       227236, 175932, 227344,
                                                                       119855, 119981, 177780,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 229828, 0, 3,
                                                                       227452, 176100, 227668,
                                                                       120233, 120443, 177948,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 230188, 0, 3,
                                                                       227668, 176268, 227884,
                                                                       120443, 120653, 178228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 230548, 0, 3,
                                                                       227884, 176436, 228100,
                                                                       120653, 120863, 178508,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 230908, 0, 3,
                                                                       228100, 176604, 228316,
                                                                       120863, 121073, 178788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 231268, 0, 3,
                                                                       228316, 176772, 228532,
                                                                       121073, 121283, 179068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 231628, 0, 3,
                                                                       228532, 176940, 228748,
                                                                       121283, 121493, 179348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 231988, 0, 3,
                                                                       228748, 177108, 228964,
                                                                       121493, 121703, 179628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 232348, 0, 3,
                                                                       228964, 177276, 229180,
                                                                       121703, 121913, 179908,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 232708, 0, 3,
                                                                       229180, 177444, 229396,
                                                                       121913, 122123, 180188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 233068, 0, 3,
                                                                       229396, 177612, 229612,
                                                                       122123, 122333, 180468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 233428, 0, 3,
                                                                       229828, 177948, 230188,
                                                                       122753, 123068, 180748,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 233968, 0, 3,
                                                                       230188, 178228, 230548,
                                                                       123068, 123383, 181168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 234508, 0, 3,
                                                                       230548, 178508, 230908,
                                                                       123383, 123698, 181588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 235048, 0, 3,
                                                                       230908, 178788, 231268,
                                                                       123698, 124013, 182008,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 235588, 0, 3,
                                                                       231268, 179068, 231628,
                                                                       124013, 124328, 182428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 236128, 0, 3,
                                                                       231628, 179348, 231988,
                                                                       124328, 124643, 182848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 236668, 0, 3,
                                                                       231988, 179628, 232348,
                                                                       124643, 124958, 183268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 237208, 0, 3,
                                                                       232348, 179908, 232708,
                                                                       124958, 125273, 183688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 237748, 0, 3,
                                                                       232708, 180188, 233068,
                                                                       125273, 125588, 184108,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 238288, 0, 3,
                                                                       233428, 180748, 233968,
                                                                       126218, 126659, 184528,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 239044, 0, 3,
                                                                       233968, 181168, 234508,
                                                                       126659, 127100, 185116,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 239800, 0, 3,
                                                                       234508, 181588, 235048,
                                                                       127100, 127541, 185704,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 240556, 0, 3,
                                                                       235048, 182008, 235588,
                                                                       127541, 127982, 186292,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 241312, 0, 3,
                                                                       235588, 182428, 236128,
                                                                       127982, 128423, 186880,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 242068, 0, 3,
                                                                       236128, 182848, 236668,
                                                                       128423, 128864, 187468,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 242824, 0, 3,
                                                                       236668, 183268, 237208,
                                                                       128864, 129305, 188056,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 243580, 0, 3,
                                                                       237208, 183688, 237748,
                                                                       129305, 129746, 188644,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 244336, 0, 3,
                                                                       238288, 184528, 239044,
                                                                       130628, 131216, 189232,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 245344, 0, 3,
                                                                       239044, 185116, 239800,
                                                                       131216, 131804, 190016,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 246352, 0, 3,
                                                                       239800, 185704, 240556,
                                                                       131804, 132392, 190800,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 247360, 0, 3,
                                                                       240556, 186292, 241312,
                                                                       132392, 132980, 191584,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 248368, 0, 3,
                                                                       241312, 186880, 242068,
                                                                       132980, 133568, 192368,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 249376, 0, 3,
                                                                       242068, 187468, 242824,
                                                                       133568, 134156, 193152,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 250384, 0, 3,
                                                                       242824, 188056, 243580,
                                                                       134156, 134744, 193936,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 251392, 0, 3,
                                                                       244336, 189232, 245344,
                                                                       135920, 136676, 194720,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 252688, 0, 3,
                                                                       245344, 190016, 246352,
                                                                       136676, 137432, 195728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 253984, 0, 3,
                                                                       246352, 190800, 247360,
                                                                       137432, 138188, 196736,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 255280, 0, 3,
                                                                       247360, 191584, 248368,
                                                                       138188, 138944, 197744,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 256576, 0, 3,
                                                                       248368, 192368, 249376,
                                                                       138944, 139700, 198752,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 257872, 0, 3,
                                                                       249376, 193152, 250384,
                                                                       139700, 140456, 199760,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 259168, 0, 3,
                                                                       251392, 194720, 252688,
                                                                       141968, 142913, 200768,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 260788, 0, 3,
                                                                       252688, 195728, 253984,
                                                                       142913, 143858, 202028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 262408, 0, 3,
                                                                       253984, 196736, 255280,
                                                                       143858, 144803, 203288,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 264028, 0, 3,
                                                                       255280, 197744, 256576,
                                                                       144803, 145748, 204548,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 265648, 0, 3,
                                                                       256576, 198752, 257872,
                                                                       145748, 146693, 205808,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 267268, 0, 3,
                                                                       259168, 200768, 260788,
                                                                       148583, 149738, 207068,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 269248, 0, 3,
                                                                       260788, 202028, 262408,
                                                                       149738, 150893, 208608,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 271228, 0, 3,
                                                                       262408, 203288, 264028,
                                                                       150893, 152048, 210148,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 273208, 0, 3,
                                                                       264028, 204548, 265648,
                                                                       152048, 153203, 211688,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 275188, 0, 3,
                                                                       267268, 207068, 269248,
                                                                       155513, 156899, 213228,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 277564, 0, 3,
                                                                       269248, 208608, 271228,
                                                                       156899, 158285, 215076,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 279940, 0, 3,
                                                                       271228, 210148, 273208,
                                                                       158285, 159671, 216924,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 282316, 0, 3,
                                                                       275188, 213228, 277564,
                                                                       162443, 164081, 218772,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 285124, 0, 3,
                                                                       277564, 215076, 279940,
                                                                       164081, 165719, 220956,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsk_three_center_electron_repulsion_0(buffer, 287932, 0, 3,
                                                                       282316, 218772, 285124,
                                                                       168995, 170906, 223140,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 291208, 244336, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 292636, 251392, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 294472, 259168, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 296767, 267268, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 299572, 275188, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 302938, 282316, 2808, ncols);

                    simdfunc::contract_primitives(buffer, 306916, 287932, 3276, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 292216, 291208, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 293932, 292636, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 296092, 294472, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 298747, 296767, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 301948, 299572, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 305746, 302938, 78, 1, nmax);

        simdtrf::transform_k_inner(buffer, 310192, 306916, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 311557, 292216, 293932, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 312817, 293932, 296092, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 314437, 296092, 298747, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 316462, 298747, 301948, 15,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 318937, 301948, 305746, 15,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 321907, 305746, 310192, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 325417, 311557, 312817, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 327937, 312817, 314437, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 331177, 314437, 316462, 15,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 335227, 316462, 318937, 15,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 340177, 318937, 321907, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 346117, 325417, 327937, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 350317, 327937, 331177, 15,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 355717, 331177, 335227, 15,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 362467, 335227, 340177, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 370717, 346117, 350317, 15,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 377017, 350317, 355717, 15,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 385117, 355717, 362467, 15,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 395242, 370717, 377017, 15,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 404062, 377017, 385117, 15,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 415402, 395242, 404062, 15,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 427162, 415402, 28, 15, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 427162, 195, nmax);
    }

    for (size_t m = 0; m < 2535; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
