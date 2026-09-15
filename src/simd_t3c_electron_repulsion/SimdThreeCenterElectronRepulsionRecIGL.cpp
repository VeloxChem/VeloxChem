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


#include "SimdThreeCenterElectronRepulsionRecIGL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
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
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_igl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_igl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 295579, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1989 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 295579, 239533, 13138, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3297, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3300, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3303, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3306, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3309, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3312, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3315, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3318, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3321, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3324, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3327, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3330, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3333, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3336, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3339, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3342, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3345, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3348, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3357, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3366, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3375, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3384, 3, 14, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3393, 3, 15, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3402, 3, 16, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3411, 3, 17, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3420, 3, 18, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3429, 3, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3438, 3, 20, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3447, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3456, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3465, 3, 23, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3474, 3, 24, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3483, 3, 25, 78,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3492, 3, 33, 93,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3510, 3, 36, 99,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3528, 3, 39, 105,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3546, 3, 42, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3564, 3, 45, 117,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3582, 3, 48, 123,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3600, 3, 51, 129,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3618, 3, 54, 135,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3636, 3, 57, 141,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3654, 3, 60, 147,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3672, 3, 63, 153,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3690, 3, 66, 159,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3708, 3, 69, 165,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3726, 3, 72, 171,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3744, 3, 75, 177,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3762, 3, 93, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3792, 3, 99, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3822, 3, 105, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3852, 3, 111, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3882, 3, 117, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3912, 3, 123, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3942, 3, 129, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3972, 3, 135, 273,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4002, 3, 141, 283,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4032, 3, 147, 293,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4062, 3, 153, 303,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4092, 3, 159, 313,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4122, 3, 165, 323,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4152, 3, 171, 333,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4182, 3, 203, 373,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4227, 3, 213, 388,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4272, 3, 223, 403,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4317, 3, 233, 418,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4362, 3, 243, 433,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4407, 3, 253, 448,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4452, 3, 263, 463,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4497, 3, 273, 478,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4542, 3, 283, 493,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4587, 3, 293, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4632, 3, 303, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4677, 3, 313, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4722, 3, 323, 553,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4767, 3, 373, 610,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4830, 3, 388, 631,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4893, 3, 403, 652,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4956, 3, 418, 673,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5019, 3, 433, 694,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5082, 3, 448, 715,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5145, 3, 463, 736,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5208, 3, 478, 757,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5271, 3, 493, 778,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5334, 3, 508, 799,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5397, 3, 523, 820,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5460, 3, 538, 841,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5523, 3, 610, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5607, 3, 631, 946,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5691, 3, 652, 974,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5775, 3, 673,
                                                                       1002, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5859, 3, 694,
                                                                       1030, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5943, 3, 715,
                                                                       1058, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6027, 3, 736,
                                                                       1086, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6111, 3, 757,
                                                                       1114, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6195, 3, 778,
                                                                       1142, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6279, 3, 799,
                                                                       1170, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6363, 3, 820,
                                                                       1198, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6447, 3, 918,
                                                                       1298, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6555, 3, 946,
                                                                       1334, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6663, 3, 974,
                                                                       1370, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6771, 3, 1002,
                                                                       1406, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6879, 3, 1030,
                                                                       1442, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6987, 3, 1058,
                                                                       1478, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7095, 3, 1086,
                                                                       1514, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7203, 3, 1114,
                                                                       1550, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7311, 3, 1142,
                                                                       1586, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7419, 3, 1170,
                                                                       1622, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7527, 3, 1298,
                                                                       1748, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7662, 3, 1334,
                                                                       1793, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7797, 3, 1370,
                                                                       1838, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7932, 3, 1406,
                                                                       1883, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8067, 3, 1442,
                                                                       1928, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8202, 3, 1478,
                                                                       1973, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8337, 3, 1514,
                                                                       2018, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8472, 3, 1550,
                                                                       2063, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8607, 3, 1586,
                                                                       2108, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8742, 3, 1748,
                                                                       2263, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8907, 3, 1793,
                                                                       2318, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9072, 3, 1838,
                                                                       2373, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9237, 3, 1883,
                                                                       2428, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9402, 3, 1928,
                                                                       2483, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9567, 3, 1973,
                                                                       2538, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9732, 3, 2018,
                                                                       2593, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9897, 3, 2063,
                                                                       2648, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10062, 3, 2263,
                                                                       2835, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10260, 3, 2318,
                                                                       2901, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10458, 3, 2373,
                                                                       2967, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10656, 3, 2428,
                                                                       3033, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10854, 3, 2483,
                                                                       3099, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11052, 3, 2538,
                                                                       3165, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11250, 3, 2593,
                                                                       3231, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11448, 3, 8, 9,
                                                                       3297, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11454, 3, 9, 10,
                                                                       3300, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11460, 3, 10, 11,
                                                                       3303, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11466, 3, 11, 12,
                                                                       3306, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11472, 3, 12, 13,
                                                                       3309, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11478, 3, 13, 14,
                                                                       3312, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11484, 3, 14, 15,
                                                                       3315, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11490, 3, 15, 16,
                                                                       3318, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11496, 3, 16, 17,
                                                                       3321, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11502, 3, 17, 18,
                                                                       3324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11508, 3, 18, 19,
                                                                       3327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11514, 3, 19, 20,
                                                                       3330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11520, 3, 20, 21,
                                                                       3333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11526, 3, 21, 22,
                                                                       3336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11532, 3, 22, 23,
                                                                       3339, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11538, 3, 23, 24,
                                                                       3342, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11544, 3, 24, 25,
                                                                       3345, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11550, 0, 3,
                                                                       11448, 3297, 11454, 3348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11568, 0, 3,
                                                                       11454, 3300, 11460, 3357,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11586, 0, 3,
                                                                       11460, 3303, 11466, 3366,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11604, 0, 3,
                                                                       11466, 3306, 11472, 3375,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11622, 0, 3,
                                                                       11472, 3309, 11478, 3384,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11640, 0, 3,
                                                                       11478, 3312, 11484, 3393,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11658, 0, 3,
                                                                       11484, 3315, 11490, 3402,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11676, 0, 3,
                                                                       11490, 3318, 11496, 3411,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11694, 0, 3,
                                                                       11496, 3321, 11502, 3420,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11712, 0, 3,
                                                                       11502, 3324, 11508, 3429,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11730, 0, 3,
                                                                       11508, 3327, 11514, 3438,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11748, 0, 3,
                                                                       11514, 3330, 11520, 3447,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11766, 0, 3,
                                                                       11520, 3333, 11526, 3456,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11784, 0, 3,
                                                                       11526, 3336, 11532, 3465,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11802, 0, 3,
                                                                       11532, 3339, 11538, 3474,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11820, 0, 3,
                                                                       11538, 3342, 11544, 3483,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11838, 0, 3,
                                                                       11550, 3348, 11568, 81,
                                                                       87, 3492, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11874, 0, 3,
                                                                       11568, 3357, 11586, 87,
                                                                       93, 3510, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11910, 0, 3,
                                                                       11586, 3366, 11604, 93,
                                                                       99, 3528, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11946, 0, 3,
                                                                       11604, 3375, 11622, 99,
                                                                       105, 3546, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11982, 0, 3,
                                                                       11622, 3384, 11640, 105,
                                                                       111, 3564, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12018, 0, 3,
                                                                       11640, 3393, 11658, 111,
                                                                       117, 3582, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12054, 0, 3,
                                                                       11658, 3402, 11676, 117,
                                                                       123, 3600, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12090, 0, 3,
                                                                       11676, 3411, 11694, 123,
                                                                       129, 3618, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12126, 0, 3,
                                                                       11694, 3420, 11712, 129,
                                                                       135, 3636, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12162, 0, 3,
                                                                       11712, 3429, 11730, 135,
                                                                       141, 3654, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12198, 0, 3,
                                                                       11730, 3438, 11748, 141,
                                                                       147, 3672, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12234, 0, 3,
                                                                       11748, 3447, 11766, 147,
                                                                       153, 3690, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12270, 0, 3,
                                                                       11766, 3456, 11784, 153,
                                                                       159, 3708, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12306, 0, 3,
                                                                       11784, 3465, 11802, 159,
                                                                       165, 3726, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12342, 0, 3,
                                                                       11802, 3474, 11820, 165,
                                                                       171, 3744, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12378, 0, 3,
                                                                       11838, 3492, 11874, 183,
                                                                       193, 3762, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12438, 0, 3,
                                                                       11874, 3510, 11910, 193,
                                                                       203, 3792, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12498, 0, 3,
                                                                       11910, 3528, 11946, 203,
                                                                       213, 3822, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12558, 0, 3,
                                                                       11946, 3546, 11982, 213,
                                                                       223, 3852, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12618, 0, 3,
                                                                       11982, 3564, 12018, 223,
                                                                       233, 3882, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12678, 0, 3,
                                                                       12018, 3582, 12054, 233,
                                                                       243, 3912, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12738, 0, 3,
                                                                       12054, 3600, 12090, 243,
                                                                       253, 3942, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12798, 0, 3,
                                                                       12090, 3618, 12126, 253,
                                                                       263, 3972, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12858, 0, 3,
                                                                       12126, 3636, 12162, 263,
                                                                       273, 4002, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12918, 0, 3,
                                                                       12162, 3654, 12198, 273,
                                                                       283, 4032, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12978, 0, 3,
                                                                       12198, 3672, 12234, 283,
                                                                       293, 4062, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13038, 0, 3,
                                                                       12234, 3690, 12270, 293,
                                                                       303, 4092, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13098, 0, 3,
                                                                       12270, 3708, 12306, 303,
                                                                       313, 4122, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13158, 0, 3,
                                                                       12306, 3726, 12342, 313,
                                                                       323, 4152, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13218, 0, 3,
                                                                       12378, 3762, 12438, 343,
                                                                       358, 4182, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13308, 0, 3,
                                                                       12438, 3792, 12498, 358,
                                                                       373, 4227, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13398, 0, 3,
                                                                       12498, 3822, 12558, 373,
                                                                       388, 4272, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13488, 0, 3,
                                                                       12558, 3852, 12618, 388,
                                                                       403, 4317, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13578, 0, 3,
                                                                       12618, 3882, 12678, 403,
                                                                       418, 4362, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13668, 0, 3,
                                                                       12678, 3912, 12738, 418,
                                                                       433, 4407, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13758, 0, 3,
                                                                       12738, 3942, 12798, 433,
                                                                       448, 4452, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13848, 0, 3,
                                                                       12798, 3972, 12858, 448,
                                                                       463, 4497, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13938, 0, 3,
                                                                       12858, 4002, 12918, 463,
                                                                       478, 4542, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14028, 0, 3,
                                                                       12918, 4032, 12978, 478,
                                                                       493, 4587, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14118, 0, 3,
                                                                       12978, 4062, 13038, 493,
                                                                       508, 4632, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14208, 0, 3,
                                                                       13038, 4092, 13098, 508,
                                                                       523, 4677, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14298, 0, 3,
                                                                       13098, 4122, 13158, 523,
                                                                       538, 4722, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14388, 0, 3,
                                                                       13218, 4182, 13308, 568,
                                                                       589, 4767, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14514, 0, 3,
                                                                       13308, 4227, 13398, 589,
                                                                       610, 4830, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14640, 0, 3,
                                                                       13398, 4272, 13488, 610,
                                                                       631, 4893, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14766, 0, 3,
                                                                       13488, 4317, 13578, 631,
                                                                       652, 4956, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14892, 0, 3,
                                                                       13578, 4362, 13668, 652,
                                                                       673, 5019, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15018, 0, 3,
                                                                       13668, 4407, 13758, 673,
                                                                       694, 5082, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15144, 0, 3,
                                                                       13758, 4452, 13848, 694,
                                                                       715, 5145, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15270, 0, 3,
                                                                       13848, 4497, 13938, 715,
                                                                       736, 5208, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15396, 0, 3,
                                                                       13938, 4542, 14028, 736,
                                                                       757, 5271, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15522, 0, 3,
                                                                       14028, 4587, 14118, 757,
                                                                       778, 5334, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15648, 0, 3,
                                                                       14118, 4632, 14208, 778,
                                                                       799, 5397, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15774, 0, 3,
                                                                       14208, 4677, 14298, 799,
                                                                       820, 5460, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15900, 0, 3,
                                                                       14388, 4767, 14514, 862,
                                                                       890, 5523, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16068, 0, 3,
                                                                       14514, 4830, 14640, 890,
                                                                       918, 5607, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16236, 0, 3,
                                                                       14640, 4893, 14766, 918,
                                                                       946, 5691, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16404, 0, 3,
                                                                       14766, 4956, 14892, 946,
                                                                       974, 5775, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16572, 0, 3,
                                                                       14892, 5019, 15018, 974,
                                                                       1002, 5859, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16740, 0, 3,
                                                                       15018, 5082, 15144, 1002,
                                                                       1030, 5943, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16908, 0, 3,
                                                                       15144, 5145, 15270, 1030,
                                                                       1058, 6027, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17076, 0, 3,
                                                                       15270, 5208, 15396, 1058,
                                                                       1086, 6111, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17244, 0, 3,
                                                                       15396, 5271, 15522, 1086,
                                                                       1114, 6195, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17412, 0, 3,
                                                                       15522, 5334, 15648, 1114,
                                                                       1142, 6279, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17580, 0, 3,
                                                                       15648, 5397, 15774, 1142,
                                                                       1170, 6363, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17748, 0, 3,
                                                                       15900, 5523, 16068, 1226,
                                                                       1262, 6447, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17964, 0, 3,
                                                                       16068, 5607, 16236, 1262,
                                                                       1298, 6555, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18180, 0, 3,
                                                                       16236, 5691, 16404, 1298,
                                                                       1334, 6663, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18396, 0, 3,
                                                                       16404, 5775, 16572, 1334,
                                                                       1370, 6771, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18612, 0, 3,
                                                                       16572, 5859, 16740, 1370,
                                                                       1406, 6879, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18828, 0, 3,
                                                                       16740, 5943, 16908, 1406,
                                                                       1442, 6987, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19044, 0, 3,
                                                                       16908, 6027, 17076, 1442,
                                                                       1478, 7095, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19260, 0, 3,
                                                                       17076, 6111, 17244, 1478,
                                                                       1514, 7203, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19476, 0, 3,
                                                                       17244, 6195, 17412, 1514,
                                                                       1550, 7311, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19692, 0, 3,
                                                                       17412, 6279, 17580, 1550,
                                                                       1586, 7419, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19908, 0, 3,
                                                                       17748, 6447, 17964, 1658,
                                                                       1703, 7527, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20178, 0, 3,
                                                                       17964, 6555, 18180, 1703,
                                                                       1748, 7662, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20448, 0, 3,
                                                                       18180, 6663, 18396, 1748,
                                                                       1793, 7797, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20718, 0, 3,
                                                                       18396, 6771, 18612, 1793,
                                                                       1838, 7932, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20988, 0, 3,
                                                                       18612, 6879, 18828, 1838,
                                                                       1883, 8067, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21258, 0, 3,
                                                                       18828, 6987, 19044, 1883,
                                                                       1928, 8202, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21528, 0, 3,
                                                                       19044, 7095, 19260, 1928,
                                                                       1973, 8337, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21798, 0, 3,
                                                                       19260, 7203, 19476, 1973,
                                                                       2018, 8472, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22068, 0, 3,
                                                                       19476, 7311, 19692, 2018,
                                                                       2063, 8607, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22338, 0, 3,
                                                                       19908, 7527, 20178, 2153,
                                                                       2208, 8742, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22668, 0, 3,
                                                                       20178, 7662, 20448, 2208,
                                                                       2263, 8907, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22998, 0, 3,
                                                                       20448, 7797, 20718, 2263,
                                                                       2318, 9072, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23328, 0, 3,
                                                                       20718, 7932, 20988, 2318,
                                                                       2373, 9237, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23658, 0, 3,
                                                                       20988, 8067, 21258, 2373,
                                                                       2428, 9402, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23988, 0, 3,
                                                                       21258, 8202, 21528, 2428,
                                                                       2483, 9567, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24318, 0, 3,
                                                                       21528, 8337, 21798, 2483,
                                                                       2538, 9732, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24648, 0, 3,
                                                                       21798, 8472, 22068, 2538,
                                                                       2593, 9897, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 24978, 0, 3,
                                                                       22338, 8742, 22668, 2703,
                                                                       2769, 10062, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 25374, 0, 3,
                                                                       22668, 8907, 22998, 2769,
                                                                       2835, 10260, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 25770, 0, 3,
                                                                       22998, 9072, 23328, 2835,
                                                                       2901, 10458, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26166, 0, 3,
                                                                       23328, 9237, 23658, 2901,
                                                                       2967, 10656, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26562, 0, 3,
                                                                       23658, 9402, 23988, 2967,
                                                                       3033, 10854, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26958, 0, 3,
                                                                       23988, 9567, 24318, 3033,
                                                                       3099, 11052, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27354, 0, 3,
                                                                       24318, 9732, 24648, 3099,
                                                                       3165, 11250, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27750, 3, 3297,
                                                                       3300, 11460, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27760, 3, 3300,
                                                                       3303, 11466, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27770, 3, 3303,
                                                                       3306, 11472, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27780, 3, 3306,
                                                                       3309, 11478, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27790, 3, 3309,
                                                                       3312, 11484, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27800, 3, 3312,
                                                                       3315, 11490, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27810, 3, 3315,
                                                                       3318, 11496, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27820, 3, 3318,
                                                                       3321, 11502, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27830, 3, 3321,
                                                                       3324, 11508, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27840, 3, 3324,
                                                                       3327, 11514, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27850, 3, 3327,
                                                                       3330, 11520, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27860, 3, 3330,
                                                                       3333, 11526, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27870, 3, 3333,
                                                                       3336, 11532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27880, 3, 3336,
                                                                       3339, 11538, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27890, 3, 3339,
                                                                       3342, 11544, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27900, 0, 3,
                                                                       27750, 11460, 27760, 3348,
                                                                       3357, 11586, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27930, 0, 3,
                                                                       27760, 11466, 27770, 3357,
                                                                       3366, 11604, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27960, 0, 3,
                                                                       27770, 11472, 27780, 3366,
                                                                       3375, 11622, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27990, 0, 3,
                                                                       27780, 11478, 27790, 3375,
                                                                       3384, 11640, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28020, 0, 3,
                                                                       27790, 11484, 27800, 3384,
                                                                       3393, 11658, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28050, 0, 3,
                                                                       27800, 11490, 27810, 3393,
                                                                       3402, 11676, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28080, 0, 3,
                                                                       27810, 11496, 27820, 3402,
                                                                       3411, 11694, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28110, 0, 3,
                                                                       27820, 11502, 27830, 3411,
                                                                       3420, 11712, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28140, 0, 3,
                                                                       27830, 11508, 27840, 3420,
                                                                       3429, 11730, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28170, 0, 3,
                                                                       27840, 11514, 27850, 3429,
                                                                       3438, 11748, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28200, 0, 3,
                                                                       27850, 11520, 27860, 3438,
                                                                       3447, 11766, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28230, 0, 3,
                                                                       27860, 11526, 27870, 3447,
                                                                       3456, 11784, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28260, 0, 3,
                                                                       27870, 11532, 27880, 3456,
                                                                       3465, 11802, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28290, 0, 3,
                                                                       27880, 11538, 27890, 3465,
                                                                       3474, 11820, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28320, 0, 3,
                                                                       27900, 11586, 27930, 3492,
                                                                       3510, 11910, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28380, 0, 3,
                                                                       27930, 11604, 27960, 3510,
                                                                       3528, 11946, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28440, 0, 3,
                                                                       27960, 11622, 27990, 3528,
                                                                       3546, 11982, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28500, 0, 3,
                                                                       27990, 11640, 28020, 3546,
                                                                       3564, 12018, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28560, 0, 3,
                                                                       28020, 11658, 28050, 3564,
                                                                       3582, 12054, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28620, 0, 3,
                                                                       28050, 11676, 28080, 3582,
                                                                       3600, 12090, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28680, 0, 3,
                                                                       28080, 11694, 28110, 3600,
                                                                       3618, 12126, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28740, 0, 3,
                                                                       28110, 11712, 28140, 3618,
                                                                       3636, 12162, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28800, 0, 3,
                                                                       28140, 11730, 28170, 3636,
                                                                       3654, 12198, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28860, 0, 3,
                                                                       28170, 11748, 28200, 3654,
                                                                       3672, 12234, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28920, 0, 3,
                                                                       28200, 11766, 28230, 3672,
                                                                       3690, 12270, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28980, 0, 3,
                                                                       28230, 11784, 28260, 3690,
                                                                       3708, 12306, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29040, 0, 3,
                                                                       28260, 11802, 28290, 3708,
                                                                       3726, 12342, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29100, 0, 3,
                                                                       28320, 11910, 28380, 3762,
                                                                       3792, 12498, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29200, 0, 3,
                                                                       28380, 11946, 28440, 3792,
                                                                       3822, 12558, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29300, 0, 3,
                                                                       28440, 11982, 28500, 3822,
                                                                       3852, 12618, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29400, 0, 3,
                                                                       28500, 12018, 28560, 3852,
                                                                       3882, 12678, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29500, 0, 3,
                                                                       28560, 12054, 28620, 3882,
                                                                       3912, 12738, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29600, 0, 3,
                                                                       28620, 12090, 28680, 3912,
                                                                       3942, 12798, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29700, 0, 3,
                                                                       28680, 12126, 28740, 3942,
                                                                       3972, 12858, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29800, 0, 3,
                                                                       28740, 12162, 28800, 3972,
                                                                       4002, 12918, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29900, 0, 3,
                                                                       28800, 12198, 28860, 4002,
                                                                       4032, 12978, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30000, 0, 3,
                                                                       28860, 12234, 28920, 4032,
                                                                       4062, 13038, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30100, 0, 3,
                                                                       28920, 12270, 28980, 4062,
                                                                       4092, 13098, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30200, 0, 3,
                                                                       28980, 12306, 29040, 4092,
                                                                       4122, 13158, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30300, 0, 3,
                                                                       29100, 12498, 29200, 4182,
                                                                       4227, 13398, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30450, 0, 3,
                                                                       29200, 12558, 29300, 4227,
                                                                       4272, 13488, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30600, 0, 3,
                                                                       29300, 12618, 29400, 4272,
                                                                       4317, 13578, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30750, 0, 3,
                                                                       29400, 12678, 29500, 4317,
                                                                       4362, 13668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30900, 0, 3,
                                                                       29500, 12738, 29600, 4362,
                                                                       4407, 13758, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31050, 0, 3,
                                                                       29600, 12798, 29700, 4407,
                                                                       4452, 13848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31200, 0, 3,
                                                                       29700, 12858, 29800, 4452,
                                                                       4497, 13938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31350, 0, 3,
                                                                       29800, 12918, 29900, 4497,
                                                                       4542, 14028, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31500, 0, 3,
                                                                       29900, 12978, 30000, 4542,
                                                                       4587, 14118, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31650, 0, 3,
                                                                       30000, 13038, 30100, 4587,
                                                                       4632, 14208, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31800, 0, 3,
                                                                       30100, 13098, 30200, 4632,
                                                                       4677, 14298, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 31950, 0, 3,
                                                                       30300, 13398, 30450, 4767,
                                                                       4830, 14640, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32160, 0, 3,
                                                                       30450, 13488, 30600, 4830,
                                                                       4893, 14766, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32370, 0, 3,
                                                                       30600, 13578, 30750, 4893,
                                                                       4956, 14892, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32580, 0, 3,
                                                                       30750, 13668, 30900, 4956,
                                                                       5019, 15018, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32790, 0, 3,
                                                                       30900, 13758, 31050, 5019,
                                                                       5082, 15144, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33000, 0, 3,
                                                                       31050, 13848, 31200, 5082,
                                                                       5145, 15270, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33210, 0, 3,
                                                                       31200, 13938, 31350, 5145,
                                                                       5208, 15396, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33420, 0, 3,
                                                                       31350, 14028, 31500, 5208,
                                                                       5271, 15522, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33630, 0, 3,
                                                                       31500, 14118, 31650, 5271,
                                                                       5334, 15648, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33840, 0, 3,
                                                                       31650, 14208, 31800, 5334,
                                                                       5397, 15774, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34050, 0, 3,
                                                                       31950, 14640, 32160, 5523,
                                                                       5607, 16236, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34330, 0, 3,
                                                                       32160, 14766, 32370, 5607,
                                                                       5691, 16404, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34610, 0, 3,
                                                                       32370, 14892, 32580, 5691,
                                                                       5775, 16572, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34890, 0, 3,
                                                                       32580, 15018, 32790, 5775,
                                                                       5859, 16740, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35170, 0, 3,
                                                                       32790, 15144, 33000, 5859,
                                                                       5943, 16908, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35450, 0, 3,
                                                                       33000, 15270, 33210, 5943,
                                                                       6027, 17076, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35730, 0, 3,
                                                                       33210, 15396, 33420, 6027,
                                                                       6111, 17244, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36010, 0, 3,
                                                                       33420, 15522, 33630, 6111,
                                                                       6195, 17412, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36290, 0, 3,
                                                                       33630, 15648, 33840, 6195,
                                                                       6279, 17580, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36570, 0, 3,
                                                                       34050, 16236, 34330, 6447,
                                                                       6555, 18180, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36930, 0, 3,
                                                                       34330, 16404, 34610, 6555,
                                                                       6663, 18396, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 37290, 0, 3,
                                                                       34610, 16572, 34890, 6663,
                                                                       6771, 18612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 37650, 0, 3,
                                                                       34890, 16740, 35170, 6771,
                                                                       6879, 18828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 38010, 0, 3,
                                                                       35170, 16908, 35450, 6879,
                                                                       6987, 19044, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 38370, 0, 3,
                                                                       35450, 17076, 35730, 6987,
                                                                       7095, 19260, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 38730, 0, 3,
                                                                       35730, 17244, 36010, 7095,
                                                                       7203, 19476, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 39090, 0, 3,
                                                                       36010, 17412, 36290, 7203,
                                                                       7311, 19692, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 39450, 0, 3,
                                                                       36570, 18180, 36930, 7527,
                                                                       7662, 20448, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 39900, 0, 3,
                                                                       36930, 18396, 37290, 7662,
                                                                       7797, 20718, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 40350, 0, 3,
                                                                       37290, 18612, 37650, 7797,
                                                                       7932, 20988, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 40800, 0, 3,
                                                                       37650, 18828, 38010, 7932,
                                                                       8067, 21258, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 41250, 0, 3,
                                                                       38010, 19044, 38370, 8067,
                                                                       8202, 21528, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 41700, 0, 3,
                                                                       38370, 19260, 38730, 8202,
                                                                       8337, 21798, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 42150, 0, 3,
                                                                       38730, 19476, 39090, 8337,
                                                                       8472, 22068, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 42600, 0, 3,
                                                                       39450, 20448, 39900, 8742,
                                                                       8907, 22998, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 43150, 0, 3,
                                                                       39900, 20718, 40350, 8907,
                                                                       9072, 23328, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 43700, 0, 3,
                                                                       40350, 20988, 40800, 9072,
                                                                       9237, 23658, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 44250, 0, 3,
                                                                       40800, 21258, 41250, 9237,
                                                                       9402, 23988, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 44800, 0, 3,
                                                                       41250, 21528, 41700, 9402,
                                                                       9567, 24318, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 45350, 0, 3,
                                                                       41700, 21798, 42150, 9567,
                                                                       9732, 24648, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 45900, 0, 3,
                                                                       42600, 22998, 43150,
                                                                       10062, 10260, 25770,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 46560, 0, 3,
                                                                       43150, 23328, 43700,
                                                                       10260, 10458, 26166,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 47220, 0, 3,
                                                                       43700, 23658, 44250,
                                                                       10458, 10656, 26562,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 47880, 0, 3,
                                                                       44250, 23988, 44800,
                                                                       10656, 10854, 26958,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 48540, 0, 3,
                                                                       44800, 24318, 45350,
                                                                       10854, 11052, 27354,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49200, 3, 11448,
                                                                       11454, 27750, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49215, 3, 11454,
                                                                       11460, 27760, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49230, 3, 11460,
                                                                       11466, 27770, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49245, 3, 11466,
                                                                       11472, 27780, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49260, 3, 11472,
                                                                       11478, 27790, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49275, 3, 11478,
                                                                       11484, 27800, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49290, 3, 11484,
                                                                       11490, 27810, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49305, 3, 11490,
                                                                       11496, 27820, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49320, 3, 11496,
                                                                       11502, 27830, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49335, 3, 11502,
                                                                       11508, 27840, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49350, 3, 11508,
                                                                       11514, 27850, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49365, 3, 11514,
                                                                       11520, 27860, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49380, 3, 11520,
                                                                       11526, 27870, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49395, 3, 11526,
                                                                       11532, 27880, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49410, 3, 11532,
                                                                       11538, 27890, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49425, 0, 3,
                                                                       49200, 27750, 49215,
                                                                       11550, 11568, 27900,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49470, 0, 3,
                                                                       49215, 27760, 49230,
                                                                       11568, 11586, 27930,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49515, 0, 3,
                                                                       49230, 27770, 49245,
                                                                       11586, 11604, 27960,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49560, 0, 3,
                                                                       49245, 27780, 49260,
                                                                       11604, 11622, 27990,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49605, 0, 3,
                                                                       49260, 27790, 49275,
                                                                       11622, 11640, 28020,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49650, 0, 3,
                                                                       49275, 27800, 49290,
                                                                       11640, 11658, 28050,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49695, 0, 3,
                                                                       49290, 27810, 49305,
                                                                       11658, 11676, 28080,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49740, 0, 3,
                                                                       49305, 27820, 49320,
                                                                       11676, 11694, 28110,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49785, 0, 3,
                                                                       49320, 27830, 49335,
                                                                       11694, 11712, 28140,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49830, 0, 3,
                                                                       49335, 27840, 49350,
                                                                       11712, 11730, 28170,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49875, 0, 3,
                                                                       49350, 27850, 49365,
                                                                       11730, 11748, 28200,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49920, 0, 3,
                                                                       49365, 27860, 49380,
                                                                       11748, 11766, 28230,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49965, 0, 3,
                                                                       49380, 27870, 49395,
                                                                       11766, 11784, 28260,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50010, 0, 3,
                                                                       49395, 27880, 49410,
                                                                       11784, 11802, 28290,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50055, 0, 3,
                                                                       49425, 27900, 49470,
                                                                       11838, 11874, 28320,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50145, 0, 3,
                                                                       49470, 27930, 49515,
                                                                       11874, 11910, 28380,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50235, 0, 3,
                                                                       49515, 27960, 49560,
                                                                       11910, 11946, 28440,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50325, 0, 3,
                                                                       49560, 27990, 49605,
                                                                       11946, 11982, 28500,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50415, 0, 3,
                                                                       49605, 28020, 49650,
                                                                       11982, 12018, 28560,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50505, 0, 3,
                                                                       49650, 28050, 49695,
                                                                       12018, 12054, 28620,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50595, 0, 3,
                                                                       49695, 28080, 49740,
                                                                       12054, 12090, 28680,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50685, 0, 3,
                                                                       49740, 28110, 49785,
                                                                       12090, 12126, 28740,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50775, 0, 3,
                                                                       49785, 28140, 49830,
                                                                       12126, 12162, 28800,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50865, 0, 3,
                                                                       49830, 28170, 49875,
                                                                       12162, 12198, 28860,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50955, 0, 3,
                                                                       49875, 28200, 49920,
                                                                       12198, 12234, 28920,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51045, 0, 3,
                                                                       49920, 28230, 49965,
                                                                       12234, 12270, 28980,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51135, 0, 3,
                                                                       49965, 28260, 50010,
                                                                       12270, 12306, 29040,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51225, 0, 3,
                                                                       50055, 28320, 50145,
                                                                       12378, 12438, 29100,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51375, 0, 3,
                                                                       50145, 28380, 50235,
                                                                       12438, 12498, 29200,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51525, 0, 3,
                                                                       50235, 28440, 50325,
                                                                       12498, 12558, 29300,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51675, 0, 3,
                                                                       50325, 28500, 50415,
                                                                       12558, 12618, 29400,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51825, 0, 3,
                                                                       50415, 28560, 50505,
                                                                       12618, 12678, 29500,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51975, 0, 3,
                                                                       50505, 28620, 50595,
                                                                       12678, 12738, 29600,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52125, 0, 3,
                                                                       50595, 28680, 50685,
                                                                       12738, 12798, 29700,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52275, 0, 3,
                                                                       50685, 28740, 50775,
                                                                       12798, 12858, 29800,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52425, 0, 3,
                                                                       50775, 28800, 50865,
                                                                       12858, 12918, 29900,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52575, 0, 3,
                                                                       50865, 28860, 50955,
                                                                       12918, 12978, 30000,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52725, 0, 3,
                                                                       50955, 28920, 51045,
                                                                       12978, 13038, 30100,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52875, 0, 3,
                                                                       51045, 28980, 51135,
                                                                       13038, 13098, 30200,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53025, 0, 3,
                                                                       51225, 29100, 51375,
                                                                       13218, 13308, 30300,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53250, 0, 3,
                                                                       51375, 29200, 51525,
                                                                       13308, 13398, 30450,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53475, 0, 3,
                                                                       51525, 29300, 51675,
                                                                       13398, 13488, 30600,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53700, 0, 3,
                                                                       51675, 29400, 51825,
                                                                       13488, 13578, 30750,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53925, 0, 3,
                                                                       51825, 29500, 51975,
                                                                       13578, 13668, 30900,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54150, 0, 3,
                                                                       51975, 29600, 52125,
                                                                       13668, 13758, 31050,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54375, 0, 3,
                                                                       52125, 29700, 52275,
                                                                       13758, 13848, 31200,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54600, 0, 3,
                                                                       52275, 29800, 52425,
                                                                       13848, 13938, 31350,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54825, 0, 3,
                                                                       52425, 29900, 52575,
                                                                       13938, 14028, 31500,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55050, 0, 3,
                                                                       52575, 30000, 52725,
                                                                       14028, 14118, 31650,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55275, 0, 3,
                                                                       52725, 30100, 52875,
                                                                       14118, 14208, 31800,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 55500, 0, 3,
                                                                       53025, 30300, 53250,
                                                                       14388, 14514, 31950,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 55815, 0, 3,
                                                                       53250, 30450, 53475,
                                                                       14514, 14640, 32160,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 56130, 0, 3,
                                                                       53475, 30600, 53700,
                                                                       14640, 14766, 32370,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 56445, 0, 3,
                                                                       53700, 30750, 53925,
                                                                       14766, 14892, 32580,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 56760, 0, 3,
                                                                       53925, 30900, 54150,
                                                                       14892, 15018, 32790,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57075, 0, 3,
                                                                       54150, 31050, 54375,
                                                                       15018, 15144, 33000,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57390, 0, 3,
                                                                       54375, 31200, 54600,
                                                                       15144, 15270, 33210,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57705, 0, 3,
                                                                       54600, 31350, 54825,
                                                                       15270, 15396, 33420,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58020, 0, 3,
                                                                       54825, 31500, 55050,
                                                                       15396, 15522, 33630,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58335, 0, 3,
                                                                       55050, 31650, 55275,
                                                                       15522, 15648, 33840,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 58650, 0, 3,
                                                                       55500, 31950, 55815,
                                                                       15900, 16068, 34050,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 59070, 0, 3,
                                                                       55815, 32160, 56130,
                                                                       16068, 16236, 34330,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 59490, 0, 3,
                                                                       56130, 32370, 56445,
                                                                       16236, 16404, 34610,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 59910, 0, 3,
                                                                       56445, 32580, 56760,
                                                                       16404, 16572, 34890,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 60330, 0, 3,
                                                                       56760, 32790, 57075,
                                                                       16572, 16740, 35170,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 60750, 0, 3,
                                                                       57075, 33000, 57390,
                                                                       16740, 16908, 35450,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61170, 0, 3,
                                                                       57390, 33210, 57705,
                                                                       16908, 17076, 35730,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61590, 0, 3,
                                                                       57705, 33420, 58020,
                                                                       17076, 17244, 36010,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 62010, 0, 3,
                                                                       58020, 33630, 58335,
                                                                       17244, 17412, 36290,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 62430, 0, 3,
                                                                       58650, 34050, 59070,
                                                                       17748, 17964, 36570,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 62970, 0, 3,
                                                                       59070, 34330, 59490,
                                                                       17964, 18180, 36930,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 63510, 0, 3,
                                                                       59490, 34610, 59910,
                                                                       18180, 18396, 37290,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 64050, 0, 3,
                                                                       59910, 34890, 60330,
                                                                       18396, 18612, 37650,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 64590, 0, 3,
                                                                       60330, 35170, 60750,
                                                                       18612, 18828, 38010,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 65130, 0, 3,
                                                                       60750, 35450, 61170,
                                                                       18828, 19044, 38370,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 65670, 0, 3,
                                                                       61170, 35730, 61590,
                                                                       19044, 19260, 38730,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 66210, 0, 3,
                                                                       61590, 36010, 62010,
                                                                       19260, 19476, 39090,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 66750, 0, 3,
                                                                       62430, 36570, 62970,
                                                                       19908, 20178, 39450,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 67425, 0, 3,
                                                                       62970, 36930, 63510,
                                                                       20178, 20448, 39900,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 68100, 0, 3,
                                                                       63510, 37290, 64050,
                                                                       20448, 20718, 40350,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 68775, 0, 3,
                                                                       64050, 37650, 64590,
                                                                       20718, 20988, 40800,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 69450, 0, 3,
                                                                       64590, 38010, 65130,
                                                                       20988, 21258, 41250,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 70125, 0, 3,
                                                                       65130, 38370, 65670,
                                                                       21258, 21528, 41700,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 70800, 0, 3,
                                                                       65670, 38730, 66210,
                                                                       21528, 21798, 42150,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 71475, 0, 3,
                                                                       66750, 39450, 67425,
                                                                       22338, 22668, 42600,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 72300, 0, 3,
                                                                       67425, 39900, 68100,
                                                                       22668, 22998, 43150,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 73125, 0, 3,
                                                                       68100, 40350, 68775,
                                                                       22998, 23328, 43700,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 73950, 0, 3,
                                                                       68775, 40800, 69450,
                                                                       23328, 23658, 44250,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 74775, 0, 3,
                                                                       69450, 41250, 70125,
                                                                       23658, 23988, 44800,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 75600, 0, 3,
                                                                       70125, 41700, 70800,
                                                                       23988, 24318, 45350,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 76425, 0, 3,
                                                                       71475, 42600, 72300,
                                                                       24978, 25374, 45900,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 77415, 0, 3,
                                                                       72300, 43150, 73125,
                                                                       25374, 25770, 46560,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 78405, 0, 3,
                                                                       73125, 43700, 73950,
                                                                       25770, 26166, 47220,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 79395, 0, 3,
                                                                       73950, 44250, 74775,
                                                                       26166, 26562, 47880,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 80385, 0, 3,
                                                                       74775, 44800, 75600,
                                                                       26562, 26958, 48540,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81375, 3, 27750,
                                                                       27760, 49230, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81396, 3, 27760,
                                                                       27770, 49245, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81417, 3, 27770,
                                                                       27780, 49260, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81438, 3, 27780,
                                                                       27790, 49275, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81459, 3, 27790,
                                                                       27800, 49290, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81480, 3, 27800,
                                                                       27810, 49305, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81501, 3, 27810,
                                                                       27820, 49320, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81522, 3, 27820,
                                                                       27830, 49335, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81543, 3, 27830,
                                                                       27840, 49350, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81564, 3, 27840,
                                                                       27850, 49365, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81585, 3, 27850,
                                                                       27860, 49380, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81606, 3, 27860,
                                                                       27870, 49395, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81627, 3, 27870,
                                                                       27880, 49410, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81648, 0, 3,
                                                                       81375, 49230, 81396,
                                                                       27900, 27930, 49515,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81711, 0, 3,
                                                                       81396, 49245, 81417,
                                                                       27930, 27960, 49560,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81774, 0, 3,
                                                                       81417, 49260, 81438,
                                                                       27960, 27990, 49605,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81837, 0, 3,
                                                                       81438, 49275, 81459,
                                                                       27990, 28020, 49650,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81900, 0, 3,
                                                                       81459, 49290, 81480,
                                                                       28020, 28050, 49695,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81963, 0, 3,
                                                                       81480, 49305, 81501,
                                                                       28050, 28080, 49740,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82026, 0, 3,
                                                                       81501, 49320, 81522,
                                                                       28080, 28110, 49785,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82089, 0, 3,
                                                                       81522, 49335, 81543,
                                                                       28110, 28140, 49830,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82152, 0, 3,
                                                                       81543, 49350, 81564,
                                                                       28140, 28170, 49875,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82215, 0, 3,
                                                                       81564, 49365, 81585,
                                                                       28170, 28200, 49920,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82278, 0, 3,
                                                                       81585, 49380, 81606,
                                                                       28200, 28230, 49965,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82341, 0, 3,
                                                                       81606, 49395, 81627,
                                                                       28230, 28260, 50010,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82404, 0, 3,
                                                                       81648, 49515, 81711,
                                                                       28320, 28380, 50235,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82530, 0, 3,
                                                                       81711, 49560, 81774,
                                                                       28380, 28440, 50325,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82656, 0, 3,
                                                                       81774, 49605, 81837,
                                                                       28440, 28500, 50415,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82782, 0, 3,
                                                                       81837, 49650, 81900,
                                                                       28500, 28560, 50505,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82908, 0, 3,
                                                                       81900, 49695, 81963,
                                                                       28560, 28620, 50595,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83034, 0, 3,
                                                                       81963, 49740, 82026,
                                                                       28620, 28680, 50685,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83160, 0, 3,
                                                                       82026, 49785, 82089,
                                                                       28680, 28740, 50775,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83286, 0, 3,
                                                                       82089, 49830, 82152,
                                                                       28740, 28800, 50865,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83412, 0, 3,
                                                                       82152, 49875, 82215,
                                                                       28800, 28860, 50955,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83538, 0, 3,
                                                                       82215, 49920, 82278,
                                                                       28860, 28920, 51045,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83664, 0, 3,
                                                                       82278, 49965, 82341,
                                                                       28920, 28980, 51135,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 83790, 0, 3,
                                                                       82404, 50235, 82530,
                                                                       29100, 29200, 51525,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84000, 0, 3,
                                                                       82530, 50325, 82656,
                                                                       29200, 29300, 51675,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84210, 0, 3,
                                                                       82656, 50415, 82782,
                                                                       29300, 29400, 51825,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84420, 0, 3,
                                                                       82782, 50505, 82908,
                                                                       29400, 29500, 51975,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84630, 0, 3,
                                                                       82908, 50595, 83034,
                                                                       29500, 29600, 52125,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84840, 0, 3,
                                                                       83034, 50685, 83160,
                                                                       29600, 29700, 52275,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 85050, 0, 3,
                                                                       83160, 50775, 83286,
                                                                       29700, 29800, 52425,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 85260, 0, 3,
                                                                       83286, 50865, 83412,
                                                                       29800, 29900, 52575,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 85470, 0, 3,
                                                                       83412, 50955, 83538,
                                                                       29900, 30000, 52725,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 85680, 0, 3,
                                                                       83538, 51045, 83664,
                                                                       30000, 30100, 52875,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 85890, 0, 3,
                                                                       83790, 51525, 84000,
                                                                       30300, 30450, 53475,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 86205, 0, 3,
                                                                       84000, 51675, 84210,
                                                                       30450, 30600, 53700,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 86520, 0, 3,
                                                                       84210, 51825, 84420,
                                                                       30600, 30750, 53925,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 86835, 0, 3,
                                                                       84420, 51975, 84630,
                                                                       30750, 30900, 54150,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 87150, 0, 3,
                                                                       84630, 52125, 84840,
                                                                       30900, 31050, 54375,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 87465, 0, 3,
                                                                       84840, 52275, 85050,
                                                                       31050, 31200, 54600,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 87780, 0, 3,
                                                                       85050, 52425, 85260,
                                                                       31200, 31350, 54825,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 88095, 0, 3,
                                                                       85260, 52575, 85470,
                                                                       31350, 31500, 55050,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 88410, 0, 3,
                                                                       85470, 52725, 85680,
                                                                       31500, 31650, 55275,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 88725, 0, 3,
                                                                       85890, 53475, 86205,
                                                                       31950, 32160, 56130,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 89166, 0, 3,
                                                                       86205, 53700, 86520,
                                                                       32160, 32370, 56445,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 89607, 0, 3,
                                                                       86520, 53925, 86835,
                                                                       32370, 32580, 56760,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 90048, 0, 3,
                                                                       86835, 54150, 87150,
                                                                       32580, 32790, 57075,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 90489, 0, 3,
                                                                       87150, 54375, 87465,
                                                                       32790, 33000, 57390,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 90930, 0, 3,
                                                                       87465, 54600, 87780,
                                                                       33000, 33210, 57705,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 91371, 0, 3,
                                                                       87780, 54825, 88095,
                                                                       33210, 33420, 58020,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 91812, 0, 3,
                                                                       88095, 55050, 88410,
                                                                       33420, 33630, 58335,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 92253, 0, 3,
                                                                       88725, 56130, 89166,
                                                                       34050, 34330, 59490,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 92841, 0, 3,
                                                                       89166, 56445, 89607,
                                                                       34330, 34610, 59910,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 93429, 0, 3,
                                                                       89607, 56760, 90048,
                                                                       34610, 34890, 60330,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 94017, 0, 3,
                                                                       90048, 57075, 90489,
                                                                       34890, 35170, 60750,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 94605, 0, 3,
                                                                       90489, 57390, 90930,
                                                                       35170, 35450, 61170,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 95193, 0, 3,
                                                                       90930, 57705, 91371,
                                                                       35450, 35730, 61590,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 95781, 0, 3,
                                                                       91371, 58020, 91812,
                                                                       35730, 36010, 62010,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 96369, 0, 3,
                                                                       92253, 59490, 92841,
                                                                       36570, 36930, 63510,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 97125, 0, 3,
                                                                       92841, 59910, 93429,
                                                                       36930, 37290, 64050,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 97881, 0, 3,
                                                                       93429, 60330, 94017,
                                                                       37290, 37650, 64590,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 98637, 0, 3,
                                                                       94017, 60750, 94605,
                                                                       37650, 38010, 65130,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 99393, 0, 3,
                                                                       94605, 61170, 95193,
                                                                       38010, 38370, 65670,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 100149, 0, 3,
                                                                       95193, 61590, 95781,
                                                                       38370, 38730, 66210,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 100905, 0, 3,
                                                                       96369, 63510, 97125,
                                                                       39450, 39900, 68100,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 101850, 0, 3,
                                                                       97125, 64050, 97881,
                                                                       39900, 40350, 68775,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 102795, 0, 3,
                                                                       97881, 64590, 98637,
                                                                       40350, 40800, 69450,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 103740, 0, 3,
                                                                       98637, 65130, 99393,
                                                                       40800, 41250, 70125,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 104685, 0, 3,
                                                                       99393, 65670, 100149,
                                                                       41250, 41700, 70800,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 105630, 0, 3,
                                                                       100905, 68100, 101850,
                                                                       42600, 43150, 73125,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 106785, 0, 3,
                                                                       101850, 68775, 102795,
                                                                       43150, 43700, 73950,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 107940, 0, 3,
                                                                       102795, 69450, 103740,
                                                                       43700, 44250, 74775,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 109095, 0, 3,
                                                                       103740, 70125, 104685,
                                                                       44250, 44800, 75600,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 110250, 0, 3,
                                                                       105630, 73125, 106785,
                                                                       45900, 46560, 78405,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 111636, 0, 3,
                                                                       106785, 73950, 107940,
                                                                       46560, 47220, 79395,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 113022, 0, 3,
                                                                       107940, 74775, 109095,
                                                                       47220, 47880, 80385,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114408, 3, 49200,
                                                                       49215, 81375, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114436, 3, 49215,
                                                                       49230, 81396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114464, 3, 49230,
                                                                       49245, 81417, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114492, 3, 49245,
                                                                       49260, 81438, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114520, 3, 49260,
                                                                       49275, 81459, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114548, 3, 49275,
                                                                       49290, 81480, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114576, 3, 49290,
                                                                       49305, 81501, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114604, 3, 49305,
                                                                       49320, 81522, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114632, 3, 49320,
                                                                       49335, 81543, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114660, 3, 49335,
                                                                       49350, 81564, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114688, 3, 49350,
                                                                       49365, 81585, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114716, 3, 49365,
                                                                       49380, 81606, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114744, 3, 49380,
                                                                       49395, 81627, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 114772, 0, 3,
                                                                       114408, 81375, 114436,
                                                                       49425, 49470, 81648,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 114856, 0, 3,
                                                                       114436, 81396, 114464,
                                                                       49470, 49515, 81711,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 114940, 0, 3,
                                                                       114464, 81417, 114492,
                                                                       49515, 49560, 81774,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115024, 0, 3,
                                                                       114492, 81438, 114520,
                                                                       49560, 49605, 81837,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115108, 0, 3,
                                                                       114520, 81459, 114548,
                                                                       49605, 49650, 81900,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115192, 0, 3,
                                                                       114548, 81480, 114576,
                                                                       49650, 49695, 81963,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115276, 0, 3,
                                                                       114576, 81501, 114604,
                                                                       49695, 49740, 82026,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115360, 0, 3,
                                                                       114604, 81522, 114632,
                                                                       49740, 49785, 82089,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115444, 0, 3,
                                                                       114632, 81543, 114660,
                                                                       49785, 49830, 82152,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115528, 0, 3,
                                                                       114660, 81564, 114688,
                                                                       49830, 49875, 82215,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115612, 0, 3,
                                                                       114688, 81585, 114716,
                                                                       49875, 49920, 82278,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115696, 0, 3,
                                                                       114716, 81606, 114744,
                                                                       49920, 49965, 82341,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 115780, 0, 3,
                                                                       114772, 81648, 114856,
                                                                       50055, 50145, 82404,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 115948, 0, 3,
                                                                       114856, 81711, 114940,
                                                                       50145, 50235, 82530,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116116, 0, 3,
                                                                       114940, 81774, 115024,
                                                                       50235, 50325, 82656,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116284, 0, 3,
                                                                       115024, 81837, 115108,
                                                                       50325, 50415, 82782,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116452, 0, 3,
                                                                       115108, 81900, 115192,
                                                                       50415, 50505, 82908,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116620, 0, 3,
                                                                       115192, 81963, 115276,
                                                                       50505, 50595, 83034,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116788, 0, 3,
                                                                       115276, 82026, 115360,
                                                                       50595, 50685, 83160,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116956, 0, 3,
                                                                       115360, 82089, 115444,
                                                                       50685, 50775, 83286,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 117124, 0, 3,
                                                                       115444, 82152, 115528,
                                                                       50775, 50865, 83412,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 117292, 0, 3,
                                                                       115528, 82215, 115612,
                                                                       50865, 50955, 83538,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 117460, 0, 3,
                                                                       115612, 82278, 115696,
                                                                       50955, 51045, 83664,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 117628, 0, 3,
                                                                       115780, 82404, 115948,
                                                                       51225, 51375, 83790,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 117908, 0, 3,
                                                                       115948, 82530, 116116,
                                                                       51375, 51525, 84000,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 118188, 0, 3,
                                                                       116116, 82656, 116284,
                                                                       51525, 51675, 84210,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 118468, 0, 3,
                                                                       116284, 82782, 116452,
                                                                       51675, 51825, 84420,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 118748, 0, 3,
                                                                       116452, 82908, 116620,
                                                                       51825, 51975, 84630,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 119028, 0, 3,
                                                                       116620, 83034, 116788,
                                                                       51975, 52125, 84840,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 119308, 0, 3,
                                                                       116788, 83160, 116956,
                                                                       52125, 52275, 85050,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 119588, 0, 3,
                                                                       116956, 83286, 117124,
                                                                       52275, 52425, 85260,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 119868, 0, 3,
                                                                       117124, 83412, 117292,
                                                                       52425, 52575, 85470,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 120148, 0, 3,
                                                                       117292, 83538, 117460,
                                                                       52575, 52725, 85680,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 120428, 0, 3,
                                                                       117628, 83790, 117908,
                                                                       53025, 53250, 85890,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 120848, 0, 3,
                                                                       117908, 84000, 118188,
                                                                       53250, 53475, 86205,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 121268, 0, 3,
                                                                       118188, 84210, 118468,
                                                                       53475, 53700, 86520,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 121688, 0, 3,
                                                                       118468, 84420, 118748,
                                                                       53700, 53925, 86835,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 122108, 0, 3,
                                                                       118748, 84630, 119028,
                                                                       53925, 54150, 87150,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 122528, 0, 3,
                                                                       119028, 84840, 119308,
                                                                       54150, 54375, 87465,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 122948, 0, 3,
                                                                       119308, 85050, 119588,
                                                                       54375, 54600, 87780,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 123368, 0, 3,
                                                                       119588, 85260, 119868,
                                                                       54600, 54825, 88095,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 123788, 0, 3,
                                                                       119868, 85470, 120148,
                                                                       54825, 55050, 88410,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 124208, 0, 3,
                                                                       120428, 85890, 120848,
                                                                       55500, 55815, 88725,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 124796, 0, 3,
                                                                       120848, 86205, 121268,
                                                                       55815, 56130, 89166,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 125384, 0, 3,
                                                                       121268, 86520, 121688,
                                                                       56130, 56445, 89607,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 125972, 0, 3,
                                                                       121688, 86835, 122108,
                                                                       56445, 56760, 90048,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 126560, 0, 3,
                                                                       122108, 87150, 122528,
                                                                       56760, 57075, 90489,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 127148, 0, 3,
                                                                       122528, 87465, 122948,
                                                                       57075, 57390, 90930,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 127736, 0, 3,
                                                                       122948, 87780, 123368,
                                                                       57390, 57705, 91371,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 128324, 0, 3,
                                                                       123368, 88095, 123788,
                                                                       57705, 58020, 91812,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 128912, 0, 3,
                                                                       124208, 88725, 124796,
                                                                       58650, 59070, 92253,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 129696, 0, 3,
                                                                       124796, 89166, 125384,
                                                                       59070, 59490, 92841,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 130480, 0, 3,
                                                                       125384, 89607, 125972,
                                                                       59490, 59910, 93429,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 131264, 0, 3,
                                                                       125972, 90048, 126560,
                                                                       59910, 60330, 94017,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 132048, 0, 3,
                                                                       126560, 90489, 127148,
                                                                       60330, 60750, 94605,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 132832, 0, 3,
                                                                       127148, 90930, 127736,
                                                                       60750, 61170, 95193,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 133616, 0, 3,
                                                                       127736, 91371, 128324,
                                                                       61170, 61590, 95781,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 134400, 0, 3,
                                                                       128912, 92253, 129696,
                                                                       62430, 62970, 96369,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 135408, 0, 3,
                                                                       129696, 92841, 130480,
                                                                       62970, 63510, 97125,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 136416, 0, 3,
                                                                       130480, 93429, 131264,
                                                                       63510, 64050, 97881,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 137424, 0, 3,
                                                                       131264, 94017, 132048,
                                                                       64050, 64590, 98637,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 138432, 0, 3,
                                                                       132048, 94605, 132832,
                                                                       64590, 65130, 99393,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 139440, 0, 3,
                                                                       132832, 95193, 133616,
                                                                       65130, 65670, 100149,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 140448, 0, 3,
                                                                       134400, 96369, 135408,
                                                                       66750, 67425, 100905,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 141708, 0, 3,
                                                                       135408, 97125, 136416,
                                                                       67425, 68100, 101850,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 142968, 0, 3,
                                                                       136416, 97881, 137424,
                                                                       68100, 68775, 102795,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 144228, 0, 3,
                                                                       137424, 98637, 138432,
                                                                       68775, 69450, 103740,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 145488, 0, 3,
                                                                       138432, 99393, 139440,
                                                                       69450, 70125, 104685,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 146748, 0, 3,
                                                                       140448, 100905, 141708,
                                                                       71475, 72300, 105630,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 148288, 0, 3,
                                                                       141708, 101850, 142968,
                                                                       72300, 73125, 106785,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 149828, 0, 3,
                                                                       142968, 102795, 144228,
                                                                       73125, 73950, 107940,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 151368, 0, 3,
                                                                       144228, 103740, 145488,
                                                                       73950, 74775, 109095,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 152908, 0, 3,
                                                                       146748, 105630, 148288,
                                                                       76425, 77415, 110250,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 154756, 0, 3,
                                                                       148288, 106785, 149828,
                                                                       77415, 78405, 111636,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 156604, 0, 3,
                                                                       149828, 107940, 151368,
                                                                       78405, 79395, 113022,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158452, 3, 81375,
                                                                       81396, 114464, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158488, 3, 81396,
                                                                       81417, 114492, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158524, 3, 81417,
                                                                       81438, 114520, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158560, 3, 81438,
                                                                       81459, 114548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158596, 3, 81459,
                                                                       81480, 114576, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158632, 3, 81480,
                                                                       81501, 114604, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158668, 3, 81501,
                                                                       81522, 114632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158704, 3, 81522,
                                                                       81543, 114660, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158740, 3, 81543,
                                                                       81564, 114688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158776, 3, 81564,
                                                                       81585, 114716, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158812, 3, 81585,
                                                                       81606, 114744, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 158848, 0, 3,
                                                                       158452, 114464, 158488,
                                                                       81648, 81711, 114940,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 158956, 0, 3,
                                                                       158488, 114492, 158524,
                                                                       81711, 81774, 115024,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159064, 0, 3,
                                                                       158524, 114520, 158560,
                                                                       81774, 81837, 115108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159172, 0, 3,
                                                                       158560, 114548, 158596,
                                                                       81837, 81900, 115192,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159280, 0, 3,
                                                                       158596, 114576, 158632,
                                                                       81900, 81963, 115276,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159388, 0, 3,
                                                                       158632, 114604, 158668,
                                                                       81963, 82026, 115360,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159496, 0, 3,
                                                                       158668, 114632, 158704,
                                                                       82026, 82089, 115444,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159604, 0, 3,
                                                                       158704, 114660, 158740,
                                                                       82089, 82152, 115528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159712, 0, 3,
                                                                       158740, 114688, 158776,
                                                                       82152, 82215, 115612,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159820, 0, 3,
                                                                       158776, 114716, 158812,
                                                                       82215, 82278, 115696,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 159928, 0, 3,
                                                                       158848, 114940, 158956,
                                                                       82404, 82530, 116116,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 160144, 0, 3,
                                                                       158956, 115024, 159064,
                                                                       82530, 82656, 116284,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 160360, 0, 3,
                                                                       159064, 115108, 159172,
                                                                       82656, 82782, 116452,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 160576, 0, 3,
                                                                       159172, 115192, 159280,
                                                                       82782, 82908, 116620,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 160792, 0, 3,
                                                                       159280, 115276, 159388,
                                                                       82908, 83034, 116788,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 161008, 0, 3,
                                                                       159388, 115360, 159496,
                                                                       83034, 83160, 116956,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 161224, 0, 3,
                                                                       159496, 115444, 159604,
                                                                       83160, 83286, 117124,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 161440, 0, 3,
                                                                       159604, 115528, 159712,
                                                                       83286, 83412, 117292,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 161656, 0, 3,
                                                                       159712, 115612, 159820,
                                                                       83412, 83538, 117460,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 161872, 0, 3,
                                                                       159928, 116116, 160144,
                                                                       83790, 84000, 118188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 162232, 0, 3,
                                                                       160144, 116284, 160360,
                                                                       84000, 84210, 118468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 162592, 0, 3,
                                                                       160360, 116452, 160576,
                                                                       84210, 84420, 118748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 162952, 0, 3,
                                                                       160576, 116620, 160792,
                                                                       84420, 84630, 119028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 163312, 0, 3,
                                                                       160792, 116788, 161008,
                                                                       84630, 84840, 119308,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 163672, 0, 3,
                                                                       161008, 116956, 161224,
                                                                       84840, 85050, 119588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 164032, 0, 3,
                                                                       161224, 117124, 161440,
                                                                       85050, 85260, 119868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 164392, 0, 3,
                                                                       161440, 117292, 161656,
                                                                       85260, 85470, 120148,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 164752, 0, 3,
                                                                       161872, 118188, 162232,
                                                                       85890, 86205, 121268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 165292, 0, 3,
                                                                       162232, 118468, 162592,
                                                                       86205, 86520, 121688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 165832, 0, 3,
                                                                       162592, 118748, 162952,
                                                                       86520, 86835, 122108,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 166372, 0, 3,
                                                                       162952, 119028, 163312,
                                                                       86835, 87150, 122528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 166912, 0, 3,
                                                                       163312, 119308, 163672,
                                                                       87150, 87465, 122948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 167452, 0, 3,
                                                                       163672, 119588, 164032,
                                                                       87465, 87780, 123368,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 167992, 0, 3,
                                                                       164032, 119868, 164392,
                                                                       87780, 88095, 123788,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 168532, 0, 3,
                                                                       164752, 121268, 165292,
                                                                       88725, 89166, 125384,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 169288, 0, 3,
                                                                       165292, 121688, 165832,
                                                                       89166, 89607, 125972,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 170044, 0, 3,
                                                                       165832, 122108, 166372,
                                                                       89607, 90048, 126560,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 170800, 0, 3,
                                                                       166372, 122528, 166912,
                                                                       90048, 90489, 127148,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 171556, 0, 3,
                                                                       166912, 122948, 167452,
                                                                       90489, 90930, 127736,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 172312, 0, 3,
                                                                       167452, 123368, 167992,
                                                                       90930, 91371, 128324,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 173068, 0, 3,
                                                                       168532, 125384, 169288,
                                                                       92253, 92841, 130480,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 174076, 0, 3,
                                                                       169288, 125972, 170044,
                                                                       92841, 93429, 131264,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 175084, 0, 3,
                                                                       170044, 126560, 170800,
                                                                       93429, 94017, 132048,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 176092, 0, 3,
                                                                       170800, 127148, 171556,
                                                                       94017, 94605, 132832,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 177100, 0, 3,
                                                                       171556, 127736, 172312,
                                                                       94605, 95193, 133616,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 178108, 0, 3,
                                                                       173068, 130480, 174076,
                                                                       96369, 97125, 136416,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 179404, 0, 3,
                                                                       174076, 131264, 175084,
                                                                       97125, 97881, 137424,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 180700, 0, 3,
                                                                       175084, 132048, 176092,
                                                                       97881, 98637, 138432,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 181996, 0, 3,
                                                                       176092, 132832, 177100,
                                                                       98637, 99393, 139440,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 183292, 0, 3,
                                                                       178108, 136416, 179404,
                                                                       100905, 101850, 142968,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 184912, 0, 3,
                                                                       179404, 137424, 180700,
                                                                       101850, 102795, 144228,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 186532, 0, 3,
                                                                       180700, 138432, 181996,
                                                                       102795, 103740, 145488,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 188152, 0, 3,
                                                                       183292, 142968, 184912,
                                                                       105630, 106785, 149828,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 190132, 0, 3,
                                                                       184912, 144228, 186532,
                                                                       106785, 107940, 151368,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 192112, 0, 3,
                                                                       188152, 149828, 190132,
                                                                       110250, 111636, 156604,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194488, 3, 114408,
                                                                       114436, 158452, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194533, 3, 114436,
                                                                       114464, 158488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194578, 3, 114464,
                                                                       114492, 158524, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194623, 3, 114492,
                                                                       114520, 158560, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194668, 3, 114520,
                                                                       114548, 158596, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194713, 3, 114548,
                                                                       114576, 158632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194758, 3, 114576,
                                                                       114604, 158668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194803, 3, 114604,
                                                                       114632, 158704, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194848, 3, 114632,
                                                                       114660, 158740, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194893, 3, 114660,
                                                                       114688, 158776, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194938, 3, 114688,
                                                                       114716, 158812, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 194983, 0, 3,
                                                                       194488, 158452, 194533,
                                                                       114772, 114856, 158848,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195118, 0, 3,
                                                                       194533, 158488, 194578,
                                                                       114856, 114940, 158956,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195253, 0, 3,
                                                                       194578, 158524, 194623,
                                                                       114940, 115024, 159064,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195388, 0, 3,
                                                                       194623, 158560, 194668,
                                                                       115024, 115108, 159172,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195523, 0, 3,
                                                                       194668, 158596, 194713,
                                                                       115108, 115192, 159280,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195658, 0, 3,
                                                                       194713, 158632, 194758,
                                                                       115192, 115276, 159388,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195793, 0, 3,
                                                                       194758, 158668, 194803,
                                                                       115276, 115360, 159496,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195928, 0, 3,
                                                                       194803, 158704, 194848,
                                                                       115360, 115444, 159604,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 196063, 0, 3,
                                                                       194848, 158740, 194893,
                                                                       115444, 115528, 159712,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 196198, 0, 3,
                                                                       194893, 158776, 194938,
                                                                       115528, 115612, 159820,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 196333, 0, 3,
                                                                       194983, 158848, 195118,
                                                                       115780, 115948, 159928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 196603, 0, 3,
                                                                       195118, 158956, 195253,
                                                                       115948, 116116, 160144,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 196873, 0, 3,
                                                                       195253, 159064, 195388,
                                                                       116116, 116284, 160360,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 197143, 0, 3,
                                                                       195388, 159172, 195523,
                                                                       116284, 116452, 160576,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 197413, 0, 3,
                                                                       195523, 159280, 195658,
                                                                       116452, 116620, 160792,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 197683, 0, 3,
                                                                       195658, 159388, 195793,
                                                                       116620, 116788, 161008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 197953, 0, 3,
                                                                       195793, 159496, 195928,
                                                                       116788, 116956, 161224,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 198223, 0, 3,
                                                                       195928, 159604, 196063,
                                                                       116956, 117124, 161440,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 198493, 0, 3,
                                                                       196063, 159712, 196198,
                                                                       117124, 117292, 161656,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 198763, 0, 3,
                                                                       196333, 159928, 196603,
                                                                       117628, 117908, 161872,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 199213, 0, 3,
                                                                       196603, 160144, 196873,
                                                                       117908, 118188, 162232,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 199663, 0, 3,
                                                                       196873, 160360, 197143,
                                                                       118188, 118468, 162592,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 200113, 0, 3,
                                                                       197143, 160576, 197413,
                                                                       118468, 118748, 162952,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 200563, 0, 3,
                                                                       197413, 160792, 197683,
                                                                       118748, 119028, 163312,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 201013, 0, 3,
                                                                       197683, 161008, 197953,
                                                                       119028, 119308, 163672,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 201463, 0, 3,
                                                                       197953, 161224, 198223,
                                                                       119308, 119588, 164032,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 201913, 0, 3,
                                                                       198223, 161440, 198493,
                                                                       119588, 119868, 164392,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 202363, 0, 3,
                                                                       198763, 161872, 199213,
                                                                       120428, 120848, 164752,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 203038, 0, 3,
                                                                       199213, 162232, 199663,
                                                                       120848, 121268, 165292,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 203713, 0, 3,
                                                                       199663, 162592, 200113,
                                                                       121268, 121688, 165832,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 204388, 0, 3,
                                                                       200113, 162952, 200563,
                                                                       121688, 122108, 166372,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 205063, 0, 3,
                                                                       200563, 163312, 201013,
                                                                       122108, 122528, 166912,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 205738, 0, 3,
                                                                       201013, 163672, 201463,
                                                                       122528, 122948, 167452,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 206413, 0, 3,
                                                                       201463, 164032, 201913,
                                                                       122948, 123368, 167992,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 207088, 0, 3,
                                                                       202363, 164752, 203038,
                                                                       124208, 124796, 168532,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 208033, 0, 3,
                                                                       203038, 165292, 203713,
                                                                       124796, 125384, 169288,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 208978, 0, 3,
                                                                       203713, 165832, 204388,
                                                                       125384, 125972, 170044,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 209923, 0, 3,
                                                                       204388, 166372, 205063,
                                                                       125972, 126560, 170800,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 210868, 0, 3,
                                                                       205063, 166912, 205738,
                                                                       126560, 127148, 171556,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 211813, 0, 3,
                                                                       205738, 167452, 206413,
                                                                       127148, 127736, 172312,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 212758, 0, 3,
                                                                       207088, 168532, 208033,
                                                                       128912, 129696, 173068,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 214018, 0, 3,
                                                                       208033, 169288, 208978,
                                                                       129696, 130480, 174076,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 215278, 0, 3,
                                                                       208978, 170044, 209923,
                                                                       130480, 131264, 175084,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 216538, 0, 3,
                                                                       209923, 170800, 210868,
                                                                       131264, 132048, 176092,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 217798, 0, 3,
                                                                       210868, 171556, 211813,
                                                                       132048, 132832, 177100,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 219058, 0, 3,
                                                                       212758, 173068, 214018,
                                                                       134400, 135408, 178108,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 220678, 0, 3,
                                                                       214018, 174076, 215278,
                                                                       135408, 136416, 179404,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 222298, 0, 3,
                                                                       215278, 175084, 216538,
                                                                       136416, 137424, 180700,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 223918, 0, 3,
                                                                       216538, 176092, 217798,
                                                                       137424, 138432, 181996,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 225538, 0, 3,
                                                                       219058, 178108, 220678,
                                                                       140448, 141708, 183292,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 227563, 0, 3,
                                                                       220678, 179404, 222298,
                                                                       141708, 142968, 184912,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 229588, 0, 3,
                                                                       222298, 180700, 223918,
                                                                       142968, 144228, 186532,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 231613, 0, 3,
                                                                       225538, 183292, 227563,
                                                                       146748, 148288, 188152,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 234088, 0, 3,
                                                                       227563, 184912, 229588,
                                                                       148288, 149828, 190132,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 236563, 0, 3,
                                                                       231613, 188152, 234088,
                                                                       152908, 154756, 192112,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 239533, 212758, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 241269, 219058, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 243501, 225538, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 246291, 231613, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 249701, 236563, 2970, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 240793, 239533, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 242889, 241269, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 245526, 243501, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 248766, 246291, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 252671, 249701, 66, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 253793, 240793, 242889, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 255221, 242889, 245526, 17,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 257057, 245526, 248766, 17,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 259352, 248766, 252671, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 262157, 253793, 255221, 17,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 265013, 255221, 257057, 17,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 268685, 257057, 259352, 17,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 273275, 262157, 265013, 17,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 278035, 265013, 268685, 17,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 284155, 273275, 278035, 17,
                                             nmax);

        simdtrf::transform_g_inner(buffer, 291295, 284155, 28, 17, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 291295, 153, nmax);
    }

    for (size_t m = 0; m < 1989; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
