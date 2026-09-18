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


#include "SimdThreeCenterElectronRepulsionRsRecHGH.hpp"

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

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_hgh_three_center_electron_repulsion(double               *values,
                                               const size_t          npairs,
                                               const size_t          natoms,
                                               const CBasisFunction &a_function,
                                               const CBasisFunction &b_function,
                                               const CBasisFunction &c_function,
                                               const CSimdMatrix    &coordinates,
                                               const CSimdMatrix    &c_coordinates,
                                               CSimdMatrix          &buffer,
                                               const double          omega,
                                               const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_hgh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 152637, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2178 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    simdfunc::compute_pair_exponents(a_function, b_function, coordinates, nmax);

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 152637, 101208, 11235, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fa = -b_exps[j] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pa(buffer, coordinates, 0, nmax, fa);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_erf_boys_function(buffer, coordinates, 6, 3, {1, 2, 3,
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                                                            14}, ncols, fj, i * nprim_b + j, fq,
                                                            omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 21, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
                                                        ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 69, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 75, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 78, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 81, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 84, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 87, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 90, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 93, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 96, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 99, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 105, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 111, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 114, 0, 3, 7, 8,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 120, 0, 3, 8, 9,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 126, 0, 3, 9, 10,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 10, 11,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 11, 12,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 144, 0, 3, 12, 13,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 150, 0, 3, 13, 14,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 156, 0, 3, 14, 15,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 15, 16,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 16, 17,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 174, 0, 3, 17, 18,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 180, 0, 3, 18, 19,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 186, 0, 3, 22, 23,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 23, 24,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 24, 25,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 204, 0, 3, 25, 26,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 210, 0, 3, 26, 27,
                                                                       87, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 216, 0, 3, 27, 28,
                                                                       90, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 28, 29,
                                                                       93, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 29, 30,
                                                                       96, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 234, 0, 3, 30, 31,
                                                                       99, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 240, 0, 3, 31, 32,
                                                                       102, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 246, 0, 3, 32, 33,
                                                                       105, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 252, 0, 3, 33, 34,
                                                                       108, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 36, 39,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 39, 42,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 42, 45,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 45, 48,
                                                                       132, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 48, 51,
                                                                       138, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 51, 54,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 54, 57,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 57, 60,
                                                                       156, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 60, 63,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 63, 66,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 66, 69,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 75, 78,
                                                                       186, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 78, 81,
                                                                       192, 198, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 81, 84,
                                                                       198, 204, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 84, 87,
                                                                       204, 210, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 87, 90,
                                                                       210, 216, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 90, 93,
                                                                       216, 222, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 93, 96,
                                                                       222, 228, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 96, 99,
                                                                       228, 234, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 99,
                                                                       102, 234, 240, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 102,
                                                                       105, 240, 246, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 105,
                                                                       108, 246, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 114,
                                                                       120, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 493, 0, 3, 120,
                                                                       126, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 126,
                                                                       132, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 523, 0, 3, 132,
                                                                       138, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 538, 0, 3, 138,
                                                                       144, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 553, 0, 3, 144,
                                                                       150, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 568, 0, 3, 150,
                                                                       156, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 583, 0, 3, 156,
                                                                       162, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 598, 0, 3, 162,
                                                                       168, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 613, 0, 3, 168,
                                                                       174, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 628, 0, 3, 186,
                                                                       192, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 643, 0, 3, 192,
                                                                       198, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 658, 0, 3, 198,
                                                                       204, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 673, 0, 3, 204,
                                                                       210, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 688, 0, 3, 210,
                                                                       216, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 703, 0, 3, 216,
                                                                       222, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 718, 0, 3, 222,
                                                                       228, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 733, 0, 3, 228,
                                                                       234, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 748, 0, 3, 234,
                                                                       240, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 763, 0, 3, 240,
                                                                       246, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 778, 0, 3, 258,
                                                                       268, 478, 493, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 799, 0, 3, 268,
                                                                       278, 493, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 820, 0, 3, 278,
                                                                       288, 508, 523, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 841, 0, 3, 288,
                                                                       298, 523, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 862, 0, 3, 298,
                                                                       308, 538, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 883, 0, 3, 308,
                                                                       318, 553, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 904, 0, 3, 318,
                                                                       328, 568, 583, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 925, 0, 3, 328,
                                                                       338, 583, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 946, 0, 3, 338,
                                                                       348, 598, 613, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 967, 0, 3, 368,
                                                                       378, 628, 643, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 988, 0, 3, 378,
                                                                       388, 643, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1009, 0, 3, 388,
                                                                       398, 658, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 398,
                                                                       408, 673, 688, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1051, 0, 3, 408,
                                                                       418, 688, 703, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 418,
                                                                       428, 703, 718, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1093, 0, 3, 428,
                                                                       438, 718, 733, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 438,
                                                                       448, 733, 748, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1135, 0, 3, 448,
                                                                       458, 748, 763, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 478,
                                                                       493, 778, 799, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 493,
                                                                       508, 799, 820, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 508,
                                                                       523, 820, 841, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 523,
                                                                       538, 841, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 538,
                                                                       553, 862, 883, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 553,
                                                                       568, 883, 904, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 568,
                                                                       583, 904, 925, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 583,
                                                                       598, 925, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 628,
                                                                       643, 967, 988, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 643,
                                                                       658, 988, 1009, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 658,
                                                                       673, 1009, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 673,
                                                                       688, 1030, 1051, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 688,
                                                                       703, 1051, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 703,
                                                                       718, 1072, 1093, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 718,
                                                                       733, 1093, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 733,
                                                                       748, 1114, 1135, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 778,
                                                                       799, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1640, 0, 3, 799,
                                                                       820, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1676, 0, 3, 820,
                                                                       841, 1212, 1240, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1712, 0, 3, 841,
                                                                       862, 1240, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1748, 0, 3, 862,
                                                                       883, 1268, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1784, 0, 3, 883,
                                                                       904, 1296, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1820, 0, 3, 904,
                                                                       925, 1324, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 967,
                                                                       988, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1892, 0, 3, 988,
                                                                       1009, 1408, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1928, 0, 3, 1009,
                                                                       1030, 1436, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1964, 0, 3, 1030,
                                                                       1051, 1464, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2000, 0, 3, 1051,
                                                                       1072, 1492, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2036, 0, 3, 1072,
                                                                       1093, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2072, 0, 3, 1093,
                                                                       1114, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 1156,
                                                                       1184, 1604, 1640, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1184,
                                                                       1212, 1640, 1676, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2198, 0, 3, 1212,
                                                                       1240, 1676, 1712, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2243, 0, 3, 1240,
                                                                       1268, 1712, 1748, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2288, 0, 3, 1268,
                                                                       1296, 1748, 1784, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2333, 0, 3, 1296,
                                                                       1324, 1784, 1820, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2378, 0, 3, 1380,
                                                                       1408, 1856, 1892, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2423, 0, 3, 1408,
                                                                       1436, 1892, 1928, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2468, 0, 3, 1436,
                                                                       1464, 1928, 1964, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2513, 0, 3, 1464,
                                                                       1492, 1964, 2000, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2558, 0, 3, 1492,
                                                                       1520, 2000, 2036, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2603, 0, 3, 1520,
                                                                       1548, 2036, 2072, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1604,
                                                                       1640, 2108, 2153, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2703, 0, 3, 1640,
                                                                       1676, 2153, 2198, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2758, 0, 3, 1676,
                                                                       1712, 2198, 2243, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2813, 0, 3, 1712,
                                                                       1748, 2243, 2288, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2868, 0, 3, 1748,
                                                                       1784, 2288, 2333, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2923, 0, 3, 1856,
                                                                       1892, 2378, 2423, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2978, 0, 3, 1892,
                                                                       1928, 2423, 2468, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3033, 0, 3, 1928,
                                                                       1964, 2468, 2513, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 1964,
                                                                       2000, 2513, 2558, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3143, 0, 3, 2000,
                                                                       2036, 2558, 2603, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3198, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3201, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3204, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3207, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3210, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3213, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3216, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3219, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3222, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3225, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3228, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3231, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3234, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3237, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3240, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3243, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3246, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3249, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3252, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3255, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3258, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3261, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3264, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3267, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3270, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3273, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3276, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3279, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3282, 3, 7, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3291, 3, 8, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3300, 3, 9, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3309, 3, 10, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3318, 3, 11, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3327, 3, 12, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3336, 3, 13, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3345, 3, 14, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3354, 3, 15, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3363, 3, 16, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3372, 3, 17, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3381, 3, 18, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3390, 3, 19, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3399, 3, 22, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3408, 3, 23, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3417, 3, 24, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3426, 3, 25, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3435, 3, 26, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3444, 3, 27, 90,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3453, 3, 28, 93,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3462, 3, 29, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3471, 3, 30, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3480, 3, 31, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3489, 3, 32, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3498, 3, 33, 108,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3507, 3, 34, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3516, 3, 36, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3534, 3, 39, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3552, 3, 42, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3570, 3, 45, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3588, 3, 48, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3606, 3, 51, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3624, 3, 54, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3642, 3, 57, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3660, 3, 60, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3678, 3, 63, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3696, 3, 66, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3714, 3, 69, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3732, 3, 75, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3750, 3, 78, 192,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3768, 3, 81, 198,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3786, 3, 84, 204,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3804, 3, 87, 210,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3822, 3, 90, 216,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3840, 3, 93, 222,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3858, 3, 96, 228,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3876, 3, 99, 234,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3894, 3, 102, 240,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3912, 3, 105, 246,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3930, 3, 108, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3948, 3, 114, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3978, 3, 120, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4008, 3, 126, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4038, 3, 132, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4068, 3, 138, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4098, 3, 144, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4128, 3, 150, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4158, 3, 156, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4188, 3, 162, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4218, 3, 168, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4248, 3, 174, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4278, 3, 186, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4308, 3, 192, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4338, 3, 198, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4368, 3, 204, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4398, 3, 210, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4428, 3, 216, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4458, 3, 222, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4488, 3, 228, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4518, 3, 234, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4548, 3, 240, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4578, 3, 246, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4608, 3, 258, 478,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4653, 3, 268, 493,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4698, 3, 278, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4743, 3, 288, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4788, 3, 298, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4833, 3, 308, 553,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4878, 3, 318, 568,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4923, 3, 328, 583,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4968, 3, 338, 598,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5013, 3, 348, 613,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5058, 3, 368, 628,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5103, 3, 378, 643,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5148, 3, 388, 658,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5193, 3, 398, 673,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5238, 3, 408, 688,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5283, 3, 418, 703,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5328, 3, 428, 718,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5373, 3, 438, 733,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5418, 3, 448, 748,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5463, 3, 458, 763,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5508, 3, 478, 778,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5571, 3, 493, 799,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5634, 3, 508, 820,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5697, 3, 523, 841,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5760, 3, 538, 862,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5823, 3, 553, 883,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5886, 3, 568, 904,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5949, 3, 583, 925,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6012, 3, 598, 946,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6075, 3, 628, 967,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6138, 3, 643, 988,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6201, 3, 658,
                                                                       1009, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6264, 3, 673,
                                                                       1030, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6327, 3, 688,
                                                                       1051, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6390, 3, 703,
                                                                       1072, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6453, 3, 718,
                                                                       1093, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6516, 3, 733,
                                                                       1114, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6579, 3, 748,
                                                                       1135, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6642, 3, 778,
                                                                       1156, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6726, 3, 799,
                                                                       1184, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6810, 3, 820,
                                                                       1212, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6894, 3, 841,
                                                                       1240, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6978, 3, 862,
                                                                       1268, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7062, 3, 883,
                                                                       1296, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7146, 3, 904,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7230, 3, 925,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7314, 3, 967,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7398, 3, 988,
                                                                       1408, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7482, 3, 1009,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7566, 3, 1030,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7650, 3, 1051,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7734, 3, 1072,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7818, 3, 1093,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7902, 3, 1114,
                                                                       1576, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7986, 3, 1156,
                                                                       1604, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8094, 3, 1184,
                                                                       1640, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8202, 3, 1212,
                                                                       1676, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8310, 3, 1240,
                                                                       1712, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8418, 3, 1268,
                                                                       1748, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8526, 3, 1296,
                                                                       1784, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8634, 3, 1324,
                                                                       1820, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8742, 3, 1380,
                                                                       1856, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8850, 3, 1408,
                                                                       1892, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8958, 3, 1436,
                                                                       1928, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9066, 3, 1464,
                                                                       1964, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9174, 3, 1492,
                                                                       2000, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9282, 3, 1520,
                                                                       2036, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9390, 3, 1548,
                                                                       2072, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9498, 3, 1604,
                                                                       2108, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9633, 3, 1640,
                                                                       2153, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9768, 3, 1676,
                                                                       2198, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9903, 3, 1712,
                                                                       2243, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10038, 3, 1748,
                                                                       2288, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10173, 3, 1784,
                                                                       2333, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10308, 3, 1856,
                                                                       2378, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10443, 3, 1892,
                                                                       2423, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10578, 3, 1928,
                                                                       2468, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10713, 3, 1964,
                                                                       2513, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10848, 3, 2000,
                                                                       2558, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10983, 3, 2036,
                                                                       2603, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11118, 3, 2108,
                                                                       2648, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11283, 3, 2153,
                                                                       2703, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11448, 3, 2198,
                                                                       2758, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11613, 3, 2243,
                                                                       2813, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11778, 3, 2288,
                                                                       2868, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11943, 3, 2378,
                                                                       2923, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12108, 3, 2423,
                                                                       2978, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12273, 3, 2468,
                                                                       3033, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12438, 3, 2513,
                                                                       3088, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12603, 3, 2558,
                                                                       3143, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12768, 3, 7, 8,
                                                                       3204, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12774, 3, 8, 9,
                                                                       3207, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12780, 3, 9, 10,
                                                                       3210, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12786, 3, 10, 11,
                                                                       3213, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12792, 3, 11, 12,
                                                                       3216, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12798, 3, 12, 13,
                                                                       3219, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12804, 3, 13, 14,
                                                                       3222, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12810, 3, 14, 15,
                                                                       3225, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12816, 3, 15, 16,
                                                                       3228, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12822, 3, 16, 17,
                                                                       3231, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12828, 3, 17, 18,
                                                                       3234, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12834, 3, 18, 19,
                                                                       3237, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12840, 3, 22, 23,
                                                                       3246, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12846, 3, 23, 24,
                                                                       3249, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12852, 3, 24, 25,
                                                                       3252, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12858, 3, 25, 26,
                                                                       3255, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12864, 3, 26, 27,
                                                                       3258, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12870, 3, 27, 28,
                                                                       3261, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12876, 3, 28, 29,
                                                                       3264, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12882, 3, 29, 30,
                                                                       3267, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12888, 3, 30, 31,
                                                                       3270, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12894, 3, 31, 32,
                                                                       3273, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12900, 3, 32, 33,
                                                                       3276, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12906, 3, 33, 34,
                                                                       3279, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12912, 0, 3,
                                                                       12768, 3204, 12774, 3300,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12930, 0, 3,
                                                                       12774, 3207, 12780, 3309,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12948, 0, 3,
                                                                       12780, 3210, 12786, 3318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12966, 0, 3,
                                                                       12786, 3213, 12792, 3327,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12984, 0, 3,
                                                                       12792, 3216, 12798, 3336,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13002, 0, 3,
                                                                       12798, 3219, 12804, 3345,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13020, 0, 3,
                                                                       12804, 3222, 12810, 3354,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13038, 0, 3,
                                                                       12810, 3225, 12816, 3363,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13056, 0, 3,
                                                                       12816, 3228, 12822, 3372,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13074, 0, 3,
                                                                       12822, 3231, 12828, 3381,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13092, 0, 3,
                                                                       12828, 3234, 12834, 3390,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13110, 0, 3,
                                                                       12840, 3246, 12846, 3417,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13128, 0, 3,
                                                                       12846, 3249, 12852, 3426,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13146, 0, 3,
                                                                       12852, 3252, 12858, 3435,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13164, 0, 3,
                                                                       12858, 3255, 12864, 3444,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13182, 0, 3,
                                                                       12864, 3258, 12870, 3453,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13200, 0, 3,
                                                                       12870, 3261, 12876, 3462,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13218, 0, 3,
                                                                       12876, 3264, 12882, 3471,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13236, 0, 3,
                                                                       12882, 3267, 12888, 3480,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13254, 0, 3,
                                                                       12888, 3270, 12894, 3489,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13272, 0, 3,
                                                                       12894, 3273, 12900, 3498,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13290, 0, 3,
                                                                       12900, 3276, 12906, 3507,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13308, 0, 3,
                                                                       12912, 3300, 12930, 114,
                                                                       120, 3552, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13344, 0, 3,
                                                                       12930, 3309, 12948, 120,
                                                                       126, 3570, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13380, 0, 3,
                                                                       12948, 3318, 12966, 126,
                                                                       132, 3588, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13416, 0, 3,
                                                                       12966, 3327, 12984, 132,
                                                                       138, 3606, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13452, 0, 3,
                                                                       12984, 3336, 13002, 138,
                                                                       144, 3624, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13488, 0, 3,
                                                                       13002, 3345, 13020, 144,
                                                                       150, 3642, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13524, 0, 3,
                                                                       13020, 3354, 13038, 150,
                                                                       156, 3660, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13560, 0, 3,
                                                                       13038, 3363, 13056, 156,
                                                                       162, 3678, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13596, 0, 3,
                                                                       13056, 3372, 13074, 162,
                                                                       168, 3696, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13632, 0, 3,
                                                                       13074, 3381, 13092, 168,
                                                                       174, 3714, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13668, 0, 3,
                                                                       13110, 3417, 13128, 186,
                                                                       192, 3768, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13704, 0, 3,
                                                                       13128, 3426, 13146, 192,
                                                                       198, 3786, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13740, 0, 3,
                                                                       13146, 3435, 13164, 198,
                                                                       204, 3804, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13776, 0, 3,
                                                                       13164, 3444, 13182, 204,
                                                                       210, 3822, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13812, 0, 3,
                                                                       13182, 3453, 13200, 210,
                                                                       216, 3840, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13848, 0, 3,
                                                                       13200, 3462, 13218, 216,
                                                                       222, 3858, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13884, 0, 3,
                                                                       13218, 3471, 13236, 222,
                                                                       228, 3876, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13920, 0, 3,
                                                                       13236, 3480, 13254, 228,
                                                                       234, 3894, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13956, 0, 3,
                                                                       13254, 3489, 13272, 234,
                                                                       240, 3912, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13992, 0, 3,
                                                                       13272, 3498, 13290, 240,
                                                                       246, 3930, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14028, 0, 3,
                                                                       13308, 3552, 13344, 258,
                                                                       268, 4008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14088, 0, 3,
                                                                       13344, 3570, 13380, 268,
                                                                       278, 4038, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14148, 0, 3,
                                                                       13380, 3588, 13416, 278,
                                                                       288, 4068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14208, 0, 3,
                                                                       13416, 3606, 13452, 288,
                                                                       298, 4098, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14268, 0, 3,
                                                                       13452, 3624, 13488, 298,
                                                                       308, 4128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14328, 0, 3,
                                                                       13488, 3642, 13524, 308,
                                                                       318, 4158, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14388, 0, 3,
                                                                       13524, 3660, 13560, 318,
                                                                       328, 4188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14448, 0, 3,
                                                                       13560, 3678, 13596, 328,
                                                                       338, 4218, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14508, 0, 3,
                                                                       13596, 3696, 13632, 338,
                                                                       348, 4248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14568, 0, 3,
                                                                       13668, 3768, 13704, 368,
                                                                       378, 4338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14628, 0, 3,
                                                                       13704, 3786, 13740, 378,
                                                                       388, 4368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14688, 0, 3,
                                                                       13740, 3804, 13776, 388,
                                                                       398, 4398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14748, 0, 3,
                                                                       13776, 3822, 13812, 398,
                                                                       408, 4428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14808, 0, 3,
                                                                       13812, 3840, 13848, 408,
                                                                       418, 4458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14868, 0, 3,
                                                                       13848, 3858, 13884, 418,
                                                                       428, 4488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14928, 0, 3,
                                                                       13884, 3876, 13920, 428,
                                                                       438, 4518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14988, 0, 3,
                                                                       13920, 3894, 13956, 438,
                                                                       448, 4548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15048, 0, 3,
                                                                       13956, 3912, 13992, 448,
                                                                       458, 4578, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15108, 0, 3,
                                                                       14028, 4008, 14088, 478,
                                                                       493, 4698, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15198, 0, 3,
                                                                       14088, 4038, 14148, 493,
                                                                       508, 4743, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15288, 0, 3,
                                                                       14148, 4068, 14208, 508,
                                                                       523, 4788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15378, 0, 3,
                                                                       14208, 4098, 14268, 523,
                                                                       538, 4833, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15468, 0, 3,
                                                                       14268, 4128, 14328, 538,
                                                                       553, 4878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15558, 0, 3,
                                                                       14328, 4158, 14388, 553,
                                                                       568, 4923, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15648, 0, 3,
                                                                       14388, 4188, 14448, 568,
                                                                       583, 4968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15738, 0, 3,
                                                                       14448, 4218, 14508, 583,
                                                                       598, 5013, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15828, 0, 3,
                                                                       14568, 4338, 14628, 628,
                                                                       643, 5148, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15918, 0, 3,
                                                                       14628, 4368, 14688, 643,
                                                                       658, 5193, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16008, 0, 3,
                                                                       14688, 4398, 14748, 658,
                                                                       673, 5238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16098, 0, 3,
                                                                       14748, 4428, 14808, 673,
                                                                       688, 5283, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16188, 0, 3,
                                                                       14808, 4458, 14868, 688,
                                                                       703, 5328, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16278, 0, 3,
                                                                       14868, 4488, 14928, 703,
                                                                       718, 5373, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16368, 0, 3,
                                                                       14928, 4518, 14988, 718,
                                                                       733, 5418, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16458, 0, 3,
                                                                       14988, 4548, 15048, 733,
                                                                       748, 5463, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16548, 0, 3,
                                                                       15108, 4698, 15198, 778,
                                                                       799, 5634, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16674, 0, 3,
                                                                       15198, 4743, 15288, 799,
                                                                       820, 5697, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16800, 0, 3,
                                                                       15288, 4788, 15378, 820,
                                                                       841, 5760, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16926, 0, 3,
                                                                       15378, 4833, 15468, 841,
                                                                       862, 5823, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17052, 0, 3,
                                                                       15468, 4878, 15558, 862,
                                                                       883, 5886, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17178, 0, 3,
                                                                       15558, 4923, 15648, 883,
                                                                       904, 5949, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17304, 0, 3,
                                                                       15648, 4968, 15738, 904,
                                                                       925, 6012, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17430, 0, 3,
                                                                       15828, 5148, 15918, 967,
                                                                       988, 6201, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17556, 0, 3,
                                                                       15918, 5193, 16008, 988,
                                                                       1009, 6264, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17682, 0, 3,
                                                                       16008, 5238, 16098, 1009,
                                                                       1030, 6327, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17808, 0, 3,
                                                                       16098, 5283, 16188, 1030,
                                                                       1051, 6390, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17934, 0, 3,
                                                                       16188, 5328, 16278, 1051,
                                                                       1072, 6453, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18060, 0, 3,
                                                                       16278, 5373, 16368, 1072,
                                                                       1093, 6516, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18186, 0, 3,
                                                                       16368, 5418, 16458, 1093,
                                                                       1114, 6579, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18312, 0, 3,
                                                                       16548, 5634, 16674, 1156,
                                                                       1184, 6810, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18480, 0, 3,
                                                                       16674, 5697, 16800, 1184,
                                                                       1212, 6894, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18648, 0, 3,
                                                                       16800, 5760, 16926, 1212,
                                                                       1240, 6978, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18816, 0, 3,
                                                                       16926, 5823, 17052, 1240,
                                                                       1268, 7062, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18984, 0, 3,
                                                                       17052, 5886, 17178, 1268,
                                                                       1296, 7146, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19152, 0, 3,
                                                                       17178, 5949, 17304, 1296,
                                                                       1324, 7230, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19320, 0, 3,
                                                                       17430, 6201, 17556, 1380,
                                                                       1408, 7482, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19488, 0, 3,
                                                                       17556, 6264, 17682, 1408,
                                                                       1436, 7566, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19656, 0, 3,
                                                                       17682, 6327, 17808, 1436,
                                                                       1464, 7650, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19824, 0, 3,
                                                                       17808, 6390, 17934, 1464,
                                                                       1492, 7734, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19992, 0, 3,
                                                                       17934, 6453, 18060, 1492,
                                                                       1520, 7818, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20160, 0, 3,
                                                                       18060, 6516, 18186, 1520,
                                                                       1548, 7902, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20328, 0, 3,
                                                                       18312, 6810, 18480, 1604,
                                                                       1640, 8202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20544, 0, 3,
                                                                       18480, 6894, 18648, 1640,
                                                                       1676, 8310, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20760, 0, 3,
                                                                       18648, 6978, 18816, 1676,
                                                                       1712, 8418, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20976, 0, 3,
                                                                       18816, 7062, 18984, 1712,
                                                                       1748, 8526, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21192, 0, 3,
                                                                       18984, 7146, 19152, 1748,
                                                                       1784, 8634, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21408, 0, 3,
                                                                       19320, 7482, 19488, 1856,
                                                                       1892, 8958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21624, 0, 3,
                                                                       19488, 7566, 19656, 1892,
                                                                       1928, 9066, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21840, 0, 3,
                                                                       19656, 7650, 19824, 1928,
                                                                       1964, 9174, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22056, 0, 3,
                                                                       19824, 7734, 19992, 1964,
                                                                       2000, 9282, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22272, 0, 3,
                                                                       19992, 7818, 20160, 2000,
                                                                       2036, 9390, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22488, 0, 3,
                                                                       20328, 8202, 20544, 2108,
                                                                       2153, 9768, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22758, 0, 3,
                                                                       20544, 8310, 20760, 2153,
                                                                       2198, 9903, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23028, 0, 3,
                                                                       20760, 8418, 20976, 2198,
                                                                       2243, 10038, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23298, 0, 3,
                                                                       20976, 8526, 21192, 2243,
                                                                       2288, 10173, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23568, 0, 3,
                                                                       21408, 8958, 21624, 2378,
                                                                       2423, 10578, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23838, 0, 3,
                                                                       21624, 9066, 21840, 2423,
                                                                       2468, 10713, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24108, 0, 3,
                                                                       21840, 9174, 22056, 2468,
                                                                       2513, 10848, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24378, 0, 3,
                                                                       22056, 9282, 22272, 2513,
                                                                       2558, 10983, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24648, 0, 3,
                                                                       22488, 9768, 22758, 2648,
                                                                       2703, 11448, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24978, 0, 3,
                                                                       22758, 9903, 23028, 2703,
                                                                       2758, 11613, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25308, 0, 3,
                                                                       23028, 10038, 23298, 2758,
                                                                       2813, 11778, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25638, 0, 3,
                                                                       23568, 10578, 23838, 2923,
                                                                       2978, 12273, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25968, 0, 3,
                                                                       23838, 10713, 24108, 2978,
                                                                       3033, 12438, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 26298, 0, 3,
                                                                       24108, 10848, 24378, 3033,
                                                                       3088, 12603, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26628, 3, 3198,
                                                                       3201, 12768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26638, 3, 3201,
                                                                       3204, 12774, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26648, 3, 3204,
                                                                       3207, 12780, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26658, 3, 3207,
                                                                       3210, 12786, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26668, 3, 3210,
                                                                       3213, 12792, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26678, 3, 3213,
                                                                       3216, 12798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26688, 3, 3216,
                                                                       3219, 12804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26698, 3, 3219,
                                                                       3222, 12810, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26708, 3, 3222,
                                                                       3225, 12816, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26718, 3, 3225,
                                                                       3228, 12822, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26728, 3, 3228,
                                                                       3231, 12828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26738, 3, 3231,
                                                                       3234, 12834, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26748, 3, 3240,
                                                                       3243, 12840, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26758, 3, 3243,
                                                                       3246, 12846, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26768, 3, 3246,
                                                                       3249, 12852, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26778, 3, 3249,
                                                                       3252, 12858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26788, 3, 3252,
                                                                       3255, 12864, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26798, 3, 3255,
                                                                       3258, 12870, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26808, 3, 3258,
                                                                       3261, 12876, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26818, 3, 3261,
                                                                       3264, 12882, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26828, 3, 3264,
                                                                       3267, 12888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26838, 3, 3267,
                                                                       3270, 12894, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26848, 3, 3270,
                                                                       3273, 12900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26858, 3, 3273,
                                                                       3276, 12906, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26868, 0, 3,
                                                                       26628, 12768, 26638, 3282,
                                                                       3291, 12912, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26898, 0, 3,
                                                                       26638, 12774, 26648, 3291,
                                                                       3300, 12930, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26928, 0, 3,
                                                                       26648, 12780, 26658, 3300,
                                                                       3309, 12948, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26958, 0, 3,
                                                                       26658, 12786, 26668, 3309,
                                                                       3318, 12966, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26988, 0, 3,
                                                                       26668, 12792, 26678, 3318,
                                                                       3327, 12984, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27018, 0, 3,
                                                                       26678, 12798, 26688, 3327,
                                                                       3336, 13002, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27048, 0, 3,
                                                                       26688, 12804, 26698, 3336,
                                                                       3345, 13020, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27078, 0, 3,
                                                                       26698, 12810, 26708, 3345,
                                                                       3354, 13038, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27108, 0, 3,
                                                                       26708, 12816, 26718, 3354,
                                                                       3363, 13056, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27138, 0, 3,
                                                                       26718, 12822, 26728, 3363,
                                                                       3372, 13074, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27168, 0, 3,
                                                                       26728, 12828, 26738, 3372,
                                                                       3381, 13092, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27198, 0, 3,
                                                                       26748, 12840, 26758, 3399,
                                                                       3408, 13110, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27228, 0, 3,
                                                                       26758, 12846, 26768, 3408,
                                                                       3417, 13128, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27258, 0, 3,
                                                                       26768, 12852, 26778, 3417,
                                                                       3426, 13146, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27288, 0, 3,
                                                                       26778, 12858, 26788, 3426,
                                                                       3435, 13164, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27318, 0, 3,
                                                                       26788, 12864, 26798, 3435,
                                                                       3444, 13182, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27348, 0, 3,
                                                                       26798, 12870, 26808, 3444,
                                                                       3453, 13200, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27378, 0, 3,
                                                                       26808, 12876, 26818, 3453,
                                                                       3462, 13218, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27408, 0, 3,
                                                                       26818, 12882, 26828, 3462,
                                                                       3471, 13236, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27438, 0, 3,
                                                                       26828, 12888, 26838, 3471,
                                                                       3480, 13254, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27468, 0, 3,
                                                                       26838, 12894, 26848, 3480,
                                                                       3489, 13272, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27498, 0, 3,
                                                                       26848, 12900, 26858, 3489,
                                                                       3498, 13290, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27528, 0, 3,
                                                                       26868, 12912, 26898, 3516,
                                                                       3534, 13308, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27588, 0, 3,
                                                                       26898, 12930, 26928, 3534,
                                                                       3552, 13344, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27648, 0, 3,
                                                                       26928, 12948, 26958, 3552,
                                                                       3570, 13380, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27708, 0, 3,
                                                                       26958, 12966, 26988, 3570,
                                                                       3588, 13416, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27768, 0, 3,
                                                                       26988, 12984, 27018, 3588,
                                                                       3606, 13452, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27828, 0, 3,
                                                                       27018, 13002, 27048, 3606,
                                                                       3624, 13488, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27888, 0, 3,
                                                                       27048, 13020, 27078, 3624,
                                                                       3642, 13524, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27948, 0, 3,
                                                                       27078, 13038, 27108, 3642,
                                                                       3660, 13560, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28008, 0, 3,
                                                                       27108, 13056, 27138, 3660,
                                                                       3678, 13596, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28068, 0, 3,
                                                                       27138, 13074, 27168, 3678,
                                                                       3696, 13632, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28128, 0, 3,
                                                                       27198, 13110, 27228, 3732,
                                                                       3750, 13668, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28188, 0, 3,
                                                                       27228, 13128, 27258, 3750,
                                                                       3768, 13704, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28248, 0, 3,
                                                                       27258, 13146, 27288, 3768,
                                                                       3786, 13740, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28308, 0, 3,
                                                                       27288, 13164, 27318, 3786,
                                                                       3804, 13776, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28368, 0, 3,
                                                                       27318, 13182, 27348, 3804,
                                                                       3822, 13812, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28428, 0, 3,
                                                                       27348, 13200, 27378, 3822,
                                                                       3840, 13848, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28488, 0, 3,
                                                                       27378, 13218, 27408, 3840,
                                                                       3858, 13884, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28548, 0, 3,
                                                                       27408, 13236, 27438, 3858,
                                                                       3876, 13920, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28608, 0, 3,
                                                                       27438, 13254, 27468, 3876,
                                                                       3894, 13956, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28668, 0, 3,
                                                                       27468, 13272, 27498, 3894,
                                                                       3912, 13992, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28728, 0, 3,
                                                                       27528, 13308, 27588, 3948,
                                                                       3978, 14028, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28828, 0, 3,
                                                                       27588, 13344, 27648, 3978,
                                                                       4008, 14088, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28928, 0, 3,
                                                                       27648, 13380, 27708, 4008,
                                                                       4038, 14148, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29028, 0, 3,
                                                                       27708, 13416, 27768, 4038,
                                                                       4068, 14208, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29128, 0, 3,
                                                                       27768, 13452, 27828, 4068,
                                                                       4098, 14268, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29228, 0, 3,
                                                                       27828, 13488, 27888, 4098,
                                                                       4128, 14328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29328, 0, 3,
                                                                       27888, 13524, 27948, 4128,
                                                                       4158, 14388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29428, 0, 3,
                                                                       27948, 13560, 28008, 4158,
                                                                       4188, 14448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29528, 0, 3,
                                                                       28008, 13596, 28068, 4188,
                                                                       4218, 14508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29628, 0, 3,
                                                                       28128, 13668, 28188, 4278,
                                                                       4308, 14568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29728, 0, 3,
                                                                       28188, 13704, 28248, 4308,
                                                                       4338, 14628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29828, 0, 3,
                                                                       28248, 13740, 28308, 4338,
                                                                       4368, 14688, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29928, 0, 3,
                                                                       28308, 13776, 28368, 4368,
                                                                       4398, 14748, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30028, 0, 3,
                                                                       28368, 13812, 28428, 4398,
                                                                       4428, 14808, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30128, 0, 3,
                                                                       28428, 13848, 28488, 4428,
                                                                       4458, 14868, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30228, 0, 3,
                                                                       28488, 13884, 28548, 4458,
                                                                       4488, 14928, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30328, 0, 3,
                                                                       28548, 13920, 28608, 4488,
                                                                       4518, 14988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30428, 0, 3,
                                                                       28608, 13956, 28668, 4518,
                                                                       4548, 15048, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30528, 0, 3,
                                                                       28728, 14028, 28828, 4608,
                                                                       4653, 15108, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30678, 0, 3,
                                                                       28828, 14088, 28928, 4653,
                                                                       4698, 15198, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30828, 0, 3,
                                                                       28928, 14148, 29028, 4698,
                                                                       4743, 15288, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30978, 0, 3,
                                                                       29028, 14208, 29128, 4743,
                                                                       4788, 15378, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31128, 0, 3,
                                                                       29128, 14268, 29228, 4788,
                                                                       4833, 15468, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31278, 0, 3,
                                                                       29228, 14328, 29328, 4833,
                                                                       4878, 15558, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31428, 0, 3,
                                                                       29328, 14388, 29428, 4878,
                                                                       4923, 15648, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31578, 0, 3,
                                                                       29428, 14448, 29528, 4923,
                                                                       4968, 15738, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31728, 0, 3,
                                                                       29628, 14568, 29728, 5058,
                                                                       5103, 15828, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31878, 0, 3,
                                                                       29728, 14628, 29828, 5103,
                                                                       5148, 15918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32028, 0, 3,
                                                                       29828, 14688, 29928, 5148,
                                                                       5193, 16008, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32178, 0, 3,
                                                                       29928, 14748, 30028, 5193,
                                                                       5238, 16098, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32328, 0, 3,
                                                                       30028, 14808, 30128, 5238,
                                                                       5283, 16188, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32478, 0, 3,
                                                                       30128, 14868, 30228, 5283,
                                                                       5328, 16278, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32628, 0, 3,
                                                                       30228, 14928, 30328, 5328,
                                                                       5373, 16368, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32778, 0, 3,
                                                                       30328, 14988, 30428, 5373,
                                                                       5418, 16458, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32928, 0, 3,
                                                                       30528, 15108, 30678, 5508,
                                                                       5571, 16548, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33138, 0, 3,
                                                                       30678, 15198, 30828, 5571,
                                                                       5634, 16674, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33348, 0, 3,
                                                                       30828, 15288, 30978, 5634,
                                                                       5697, 16800, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33558, 0, 3,
                                                                       30978, 15378, 31128, 5697,
                                                                       5760, 16926, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33768, 0, 3,
                                                                       31128, 15468, 31278, 5760,
                                                                       5823, 17052, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33978, 0, 3,
                                                                       31278, 15558, 31428, 5823,
                                                                       5886, 17178, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34188, 0, 3,
                                                                       31428, 15648, 31578, 5886,
                                                                       5949, 17304, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34398, 0, 3,
                                                                       31728, 15828, 31878, 6075,
                                                                       6138, 17430, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34608, 0, 3,
                                                                       31878, 15918, 32028, 6138,
                                                                       6201, 17556, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34818, 0, 3,
                                                                       32028, 16008, 32178, 6201,
                                                                       6264, 17682, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35028, 0, 3,
                                                                       32178, 16098, 32328, 6264,
                                                                       6327, 17808, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35238, 0, 3,
                                                                       32328, 16188, 32478, 6327,
                                                                       6390, 17934, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35448, 0, 3,
                                                                       32478, 16278, 32628, 6390,
                                                                       6453, 18060, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35658, 0, 3,
                                                                       32628, 16368, 32778, 6453,
                                                                       6516, 18186, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35868, 0, 3,
                                                                       32928, 16548, 33138, 6642,
                                                                       6726, 18312, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36148, 0, 3,
                                                                       33138, 16674, 33348, 6726,
                                                                       6810, 18480, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36428, 0, 3,
                                                                       33348, 16800, 33558, 6810,
                                                                       6894, 18648, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36708, 0, 3,
                                                                       33558, 16926, 33768, 6894,
                                                                       6978, 18816, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36988, 0, 3,
                                                                       33768, 17052, 33978, 6978,
                                                                       7062, 18984, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37268, 0, 3,
                                                                       33978, 17178, 34188, 7062,
                                                                       7146, 19152, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37548, 0, 3,
                                                                       34398, 17430, 34608, 7314,
                                                                       7398, 19320, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37828, 0, 3,
                                                                       34608, 17556, 34818, 7398,
                                                                       7482, 19488, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38108, 0, 3,
                                                                       34818, 17682, 35028, 7482,
                                                                       7566, 19656, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38388, 0, 3,
                                                                       35028, 17808, 35238, 7566,
                                                                       7650, 19824, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38668, 0, 3,
                                                                       35238, 17934, 35448, 7650,
                                                                       7734, 19992, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38948, 0, 3,
                                                                       35448, 18060, 35658, 7734,
                                                                       7818, 20160, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 39228, 0, 3,
                                                                       35868, 18312, 36148, 7986,
                                                                       8094, 20328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 39588, 0, 3,
                                                                       36148, 18480, 36428, 8094,
                                                                       8202, 20544, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 39948, 0, 3,
                                                                       36428, 18648, 36708, 8202,
                                                                       8310, 20760, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40308, 0, 3,
                                                                       36708, 18816, 36988, 8310,
                                                                       8418, 20976, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40668, 0, 3,
                                                                       36988, 18984, 37268, 8418,
                                                                       8526, 21192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41028, 0, 3,
                                                                       37548, 19320, 37828, 8742,
                                                                       8850, 21408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41388, 0, 3,
                                                                       37828, 19488, 38108, 8850,
                                                                       8958, 21624, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41748, 0, 3,
                                                                       38108, 19656, 38388, 8958,
                                                                       9066, 21840, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42108, 0, 3,
                                                                       38388, 19824, 38668, 9066,
                                                                       9174, 22056, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42468, 0, 3,
                                                                       38668, 19992, 38948, 9174,
                                                                       9282, 22272, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 42828, 0, 3,
                                                                       39228, 20328, 39588, 9498,
                                                                       9633, 22488, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 43278, 0, 3,
                                                                       39588, 20544, 39948, 9633,
                                                                       9768, 22758, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 43728, 0, 3,
                                                                       39948, 20760, 40308, 9768,
                                                                       9903, 23028, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44178, 0, 3,
                                                                       40308, 20976, 40668, 9903,
                                                                       10038, 23298, ncols,
                                                                       gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44628, 0, 3,
                                                                       41028, 21408, 41388,
                                                                       10308, 10443, 23568,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 45078, 0, 3,
                                                                       41388, 21624, 41748,
                                                                       10443, 10578, 23838,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 45528, 0, 3,
                                                                       41748, 21840, 42108,
                                                                       10578, 10713, 24108,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 45978, 0, 3,
                                                                       42108, 22056, 42468,
                                                                       10713, 10848, 24378,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 46428, 0, 3,
                                                                       42828, 22488, 43278,
                                                                       11118, 11283, 24648,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 46978, 0, 3,
                                                                       43278, 22758, 43728,
                                                                       11283, 11448, 24978,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 47528, 0, 3,
                                                                       43728, 23028, 44178,
                                                                       11448, 11613, 25308,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 48078, 0, 3,
                                                                       44628, 23568, 45078,
                                                                       11943, 12108, 25638,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 48628, 0, 3,
                                                                       45078, 23838, 45528,
                                                                       12108, 12273, 25968,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 49178, 0, 3,
                                                                       45528, 24108, 45978,
                                                                       12273, 12438, 26298,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49728, 3, 12768,
                                                                       12774, 26648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49743, 3, 12774,
                                                                       12780, 26658, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49758, 3, 12780,
                                                                       12786, 26668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49773, 3, 12786,
                                                                       12792, 26678, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49788, 3, 12792,
                                                                       12798, 26688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49803, 3, 12798,
                                                                       12804, 26698, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49818, 3, 12804,
                                                                       12810, 26708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49833, 3, 12810,
                                                                       12816, 26718, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49848, 3, 12816,
                                                                       12822, 26728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49863, 3, 12822,
                                                                       12828, 26738, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49878, 3, 12840,
                                                                       12846, 26768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49893, 3, 12846,
                                                                       12852, 26778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49908, 3, 12852,
                                                                       12858, 26788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49923, 3, 12858,
                                                                       12864, 26798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49938, 3, 12864,
                                                                       12870, 26808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49953, 3, 12870,
                                                                       12876, 26818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49968, 3, 12876,
                                                                       12882, 26828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49983, 3, 12882,
                                                                       12888, 26838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49998, 3, 12888,
                                                                       12894, 26848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 50013, 3, 12894,
                                                                       12900, 26858, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50028, 0, 3,
                                                                       49728, 26648, 49743,
                                                                       12912, 12930, 26928,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50073, 0, 3,
                                                                       49743, 26658, 49758,
                                                                       12930, 12948, 26958,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50118, 0, 3,
                                                                       49758, 26668, 49773,
                                                                       12948, 12966, 26988,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50163, 0, 3,
                                                                       49773, 26678, 49788,
                                                                       12966, 12984, 27018,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50208, 0, 3,
                                                                       49788, 26688, 49803,
                                                                       12984, 13002, 27048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50253, 0, 3,
                                                                       49803, 26698, 49818,
                                                                       13002, 13020, 27078,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50298, 0, 3,
                                                                       49818, 26708, 49833,
                                                                       13020, 13038, 27108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50343, 0, 3,
                                                                       49833, 26718, 49848,
                                                                       13038, 13056, 27138,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50388, 0, 3,
                                                                       49848, 26728, 49863,
                                                                       13056, 13074, 27168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50433, 0, 3,
                                                                       49878, 26768, 49893,
                                                                       13110, 13128, 27258,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50478, 0, 3,
                                                                       49893, 26778, 49908,
                                                                       13128, 13146, 27288,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50523, 0, 3,
                                                                       49908, 26788, 49923,
                                                                       13146, 13164, 27318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50568, 0, 3,
                                                                       49923, 26798, 49938,
                                                                       13164, 13182, 27348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50613, 0, 3,
                                                                       49938, 26808, 49953,
                                                                       13182, 13200, 27378,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50658, 0, 3,
                                                                       49953, 26818, 49968,
                                                                       13200, 13218, 27408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50703, 0, 3,
                                                                       49968, 26828, 49983,
                                                                       13218, 13236, 27438,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50748, 0, 3,
                                                                       49983, 26838, 49998,
                                                                       13236, 13254, 27468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50793, 0, 3,
                                                                       49998, 26848, 50013,
                                                                       13254, 13272, 27498,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50838, 0, 3,
                                                                       50028, 26928, 50073,
                                                                       13308, 13344, 27648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50928, 0, 3,
                                                                       50073, 26958, 50118,
                                                                       13344, 13380, 27708,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51018, 0, 3,
                                                                       50118, 26988, 50163,
                                                                       13380, 13416, 27768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51108, 0, 3,
                                                                       50163, 27018, 50208,
                                                                       13416, 13452, 27828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51198, 0, 3,
                                                                       50208, 27048, 50253,
                                                                       13452, 13488, 27888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51288, 0, 3,
                                                                       50253, 27078, 50298,
                                                                       13488, 13524, 27948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51378, 0, 3,
                                                                       50298, 27108, 50343,
                                                                       13524, 13560, 28008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51468, 0, 3,
                                                                       50343, 27138, 50388,
                                                                       13560, 13596, 28068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51558, 0, 3,
                                                                       50433, 27258, 50478,
                                                                       13668, 13704, 28248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51648, 0, 3,
                                                                       50478, 27288, 50523,
                                                                       13704, 13740, 28308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51738, 0, 3,
                                                                       50523, 27318, 50568,
                                                                       13740, 13776, 28368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51828, 0, 3,
                                                                       50568, 27348, 50613,
                                                                       13776, 13812, 28428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51918, 0, 3,
                                                                       50613, 27378, 50658,
                                                                       13812, 13848, 28488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 52008, 0, 3,
                                                                       50658, 27408, 50703,
                                                                       13848, 13884, 28548,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 52098, 0, 3,
                                                                       50703, 27438, 50748,
                                                                       13884, 13920, 28608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 52188, 0, 3,
                                                                       50748, 27468, 50793,
                                                                       13920, 13956, 28668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52278, 0, 3,
                                                                       50838, 27648, 50928,
                                                                       14028, 14088, 28928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52428, 0, 3,
                                                                       50928, 27708, 51018,
                                                                       14088, 14148, 29028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52578, 0, 3,
                                                                       51018, 27768, 51108,
                                                                       14148, 14208, 29128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52728, 0, 3,
                                                                       51108, 27828, 51198,
                                                                       14208, 14268, 29228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52878, 0, 3,
                                                                       51198, 27888, 51288,
                                                                       14268, 14328, 29328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53028, 0, 3,
                                                                       51288, 27948, 51378,
                                                                       14328, 14388, 29428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53178, 0, 3,
                                                                       51378, 28008, 51468,
                                                                       14388, 14448, 29528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53328, 0, 3,
                                                                       51558, 28248, 51648,
                                                                       14568, 14628, 29828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53478, 0, 3,
                                                                       51648, 28308, 51738,
                                                                       14628, 14688, 29928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53628, 0, 3,
                                                                       51738, 28368, 51828,
                                                                       14688, 14748, 30028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53778, 0, 3,
                                                                       51828, 28428, 51918,
                                                                       14748, 14808, 30128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53928, 0, 3,
                                                                       51918, 28488, 52008,
                                                                       14808, 14868, 30228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 54078, 0, 3,
                                                                       52008, 28548, 52098,
                                                                       14868, 14928, 30328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 54228, 0, 3,
                                                                       52098, 28608, 52188,
                                                                       14928, 14988, 30428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54378, 0, 3,
                                                                       52278, 28928, 52428,
                                                                       15108, 15198, 30828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54603, 0, 3,
                                                                       52428, 29028, 52578,
                                                                       15198, 15288, 30978,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54828, 0, 3,
                                                                       52578, 29128, 52728,
                                                                       15288, 15378, 31128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55053, 0, 3,
                                                                       52728, 29228, 52878,
                                                                       15378, 15468, 31278,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55278, 0, 3,
                                                                       52878, 29328, 53028,
                                                                       15468, 15558, 31428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55503, 0, 3,
                                                                       53028, 29428, 53178,
                                                                       15558, 15648, 31578,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55728, 0, 3,
                                                                       53328, 29828, 53478,
                                                                       15828, 15918, 32028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55953, 0, 3,
                                                                       53478, 29928, 53628,
                                                                       15918, 16008, 32178,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 56178, 0, 3,
                                                                       53628, 30028, 53778,
                                                                       16008, 16098, 32328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 56403, 0, 3,
                                                                       53778, 30128, 53928,
                                                                       16098, 16188, 32478,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 56628, 0, 3,
                                                                       53928, 30228, 54078,
                                                                       16188, 16278, 32628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 56853, 0, 3,
                                                                       54078, 30328, 54228,
                                                                       16278, 16368, 32778,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57078, 0, 3,
                                                                       54378, 30828, 54603,
                                                                       16548, 16674, 33348,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57393, 0, 3,
                                                                       54603, 30978, 54828,
                                                                       16674, 16800, 33558,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57708, 0, 3,
                                                                       54828, 31128, 55053,
                                                                       16800, 16926, 33768,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58023, 0, 3,
                                                                       55053, 31278, 55278,
                                                                       16926, 17052, 33978,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58338, 0, 3,
                                                                       55278, 31428, 55503,
                                                                       17052, 17178, 34188,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58653, 0, 3,
                                                                       55728, 32028, 55953,
                                                                       17430, 17556, 34818,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58968, 0, 3,
                                                                       55953, 32178, 56178,
                                                                       17556, 17682, 35028,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 59283, 0, 3,
                                                                       56178, 32328, 56403,
                                                                       17682, 17808, 35238,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 59598, 0, 3,
                                                                       56403, 32478, 56628,
                                                                       17808, 17934, 35448,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 59913, 0, 3,
                                                                       56628, 32628, 56853,
                                                                       17934, 18060, 35658,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 60228, 0, 3,
                                                                       57078, 33348, 57393,
                                                                       18312, 18480, 36428,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 60648, 0, 3,
                                                                       57393, 33558, 57708,
                                                                       18480, 18648, 36708,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61068, 0, 3,
                                                                       57708, 33768, 58023,
                                                                       18648, 18816, 36988,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61488, 0, 3,
                                                                       58023, 33978, 58338,
                                                                       18816, 18984, 37268,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61908, 0, 3,
                                                                       58653, 34818, 58968,
                                                                       19320, 19488, 38108,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 62328, 0, 3,
                                                                       58968, 35028, 59283,
                                                                       19488, 19656, 38388,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 62748, 0, 3,
                                                                       59283, 35238, 59598,
                                                                       19656, 19824, 38668,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 63168, 0, 3,
                                                                       59598, 35448, 59913,
                                                                       19824, 19992, 38948,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 63588, 0, 3,
                                                                       60228, 36428, 60648,
                                                                       20328, 20544, 39948,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 64128, 0, 3,
                                                                       60648, 36708, 61068,
                                                                       20544, 20760, 40308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 64668, 0, 3,
                                                                       61068, 36988, 61488,
                                                                       20760, 20976, 40668,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 65208, 0, 3,
                                                                       61908, 38108, 62328,
                                                                       21408, 21624, 41748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 65748, 0, 3,
                                                                       62328, 38388, 62748,
                                                                       21624, 21840, 42108,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 66288, 0, 3,
                                                                       62748, 38668, 63168,
                                                                       21840, 22056, 42468,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 66828, 0, 3,
                                                                       63588, 39948, 64128,
                                                                       22488, 22758, 43728,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 67503, 0, 3,
                                                                       64128, 40308, 64668,
                                                                       22758, 23028, 44178,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 68178, 0, 3,
                                                                       65208, 41748, 65748,
                                                                       23568, 23838, 45528,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 68853, 0, 3,
                                                                       65748, 42108, 66288,
                                                                       23838, 24108, 45978,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 69528, 0, 3,
                                                                       66828, 43728, 67503,
                                                                       24648, 24978, 47528,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 70353, 0, 3,
                                                                       68178, 45528, 68853,
                                                                       25638, 25968, 49178,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71178, 3, 26628,
                                                                       26638, 49728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71199, 3, 26638,
                                                                       26648, 49743, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71220, 3, 26648,
                                                                       26658, 49758, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71241, 3, 26658,
                                                                       26668, 49773, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71262, 3, 26668,
                                                                       26678, 49788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71283, 3, 26678,
                                                                       26688, 49803, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71304, 3, 26688,
                                                                       26698, 49818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71325, 3, 26698,
                                                                       26708, 49833, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71346, 3, 26708,
                                                                       26718, 49848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71367, 3, 26718,
                                                                       26728, 49863, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71388, 3, 26748,
                                                                       26758, 49878, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71409, 3, 26758,
                                                                       26768, 49893, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71430, 3, 26768,
                                                                       26778, 49908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71451, 3, 26778,
                                                                       26788, 49923, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71472, 3, 26788,
                                                                       26798, 49938, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71493, 3, 26798,
                                                                       26808, 49953, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71514, 3, 26808,
                                                                       26818, 49968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71535, 3, 26818,
                                                                       26828, 49983, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71556, 3, 26828,
                                                                       26838, 49998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71577, 3, 26838,
                                                                       26848, 50013, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 71598, 0, 3,
                                                                       71178, 49728, 71199,
                                                                       26868, 26898, 50028,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 71661, 0, 3,
                                                                       71199, 49743, 71220,
                                                                       26898, 26928, 50073,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 71724, 0, 3,
                                                                       71220, 49758, 71241,
                                                                       26928, 26958, 50118,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 71787, 0, 3,
                                                                       71241, 49773, 71262,
                                                                       26958, 26988, 50163,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 71850, 0, 3,
                                                                       71262, 49788, 71283,
                                                                       26988, 27018, 50208,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 71913, 0, 3,
                                                                       71283, 49803, 71304,
                                                                       27018, 27048, 50253,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 71976, 0, 3,
                                                                       71304, 49818, 71325,
                                                                       27048, 27078, 50298,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 72039, 0, 3,
                                                                       71325, 49833, 71346,
                                                                       27078, 27108, 50343,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 72102, 0, 3,
                                                                       71346, 49848, 71367,
                                                                       27108, 27138, 50388,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 72165, 0, 3,
                                                                       71388, 49878, 71409,
                                                                       27198, 27228, 50433,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 72228, 0, 3,
                                                                       71409, 49893, 71430,
                                                                       27228, 27258, 50478,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 72291, 0, 3,
                                                                       71430, 49908, 71451,
                                                                       27258, 27288, 50523,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 72354, 0, 3,
                                                                       71451, 49923, 71472,
                                                                       27288, 27318, 50568,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 72417, 0, 3,
                                                                       71472, 49938, 71493,
                                                                       27318, 27348, 50613,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 72480, 0, 3,
                                                                       71493, 49953, 71514,
                                                                       27348, 27378, 50658,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 72543, 0, 3,
                                                                       71514, 49968, 71535,
                                                                       27378, 27408, 50703,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 72606, 0, 3,
                                                                       71535, 49983, 71556,
                                                                       27408, 27438, 50748,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 72669, 0, 3,
                                                                       71556, 49998, 71577,
                                                                       27438, 27468, 50793,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 72732, 0, 3,
                                                                       71598, 50028, 71661,
                                                                       27528, 27588, 50838,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 72858, 0, 3,
                                                                       71661, 50073, 71724,
                                                                       27588, 27648, 50928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 72984, 0, 3,
                                                                       71724, 50118, 71787,
                                                                       27648, 27708, 51018,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 73110, 0, 3,
                                                                       71787, 50163, 71850,
                                                                       27708, 27768, 51108,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 73236, 0, 3,
                                                                       71850, 50208, 71913,
                                                                       27768, 27828, 51198,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 73362, 0, 3,
                                                                       71913, 50253, 71976,
                                                                       27828, 27888, 51288,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 73488, 0, 3,
                                                                       71976, 50298, 72039,
                                                                       27888, 27948, 51378,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 73614, 0, 3,
                                                                       72039, 50343, 72102,
                                                                       27948, 28008, 51468,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 73740, 0, 3,
                                                                       72165, 50433, 72228,
                                                                       28128, 28188, 51558,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 73866, 0, 3,
                                                                       72228, 50478, 72291,
                                                                       28188, 28248, 51648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 73992, 0, 3,
                                                                       72291, 50523, 72354,
                                                                       28248, 28308, 51738,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 74118, 0, 3,
                                                                       72354, 50568, 72417,
                                                                       28308, 28368, 51828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 74244, 0, 3,
                                                                       72417, 50613, 72480,
                                                                       28368, 28428, 51918,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 74370, 0, 3,
                                                                       72480, 50658, 72543,
                                                                       28428, 28488, 52008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 74496, 0, 3,
                                                                       72543, 50703, 72606,
                                                                       28488, 28548, 52098,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 74622, 0, 3,
                                                                       72606, 50748, 72669,
                                                                       28548, 28608, 52188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 74748, 0, 3,
                                                                       72732, 50838, 72858,
                                                                       28728, 28828, 52278,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 74958, 0, 3,
                                                                       72858, 50928, 72984,
                                                                       28828, 28928, 52428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 75168, 0, 3,
                                                                       72984, 51018, 73110,
                                                                       28928, 29028, 52578,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 75378, 0, 3,
                                                                       73110, 51108, 73236,
                                                                       29028, 29128, 52728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 75588, 0, 3,
                                                                       73236, 51198, 73362,
                                                                       29128, 29228, 52878,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 75798, 0, 3,
                                                                       73362, 51288, 73488,
                                                                       29228, 29328, 53028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 76008, 0, 3,
                                                                       73488, 51378, 73614,
                                                                       29328, 29428, 53178,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 76218, 0, 3,
                                                                       73740, 51558, 73866,
                                                                       29628, 29728, 53328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 76428, 0, 3,
                                                                       73866, 51648, 73992,
                                                                       29728, 29828, 53478,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 76638, 0, 3,
                                                                       73992, 51738, 74118,
                                                                       29828, 29928, 53628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 76848, 0, 3,
                                                                       74118, 51828, 74244,
                                                                       29928, 30028, 53778,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 77058, 0, 3,
                                                                       74244, 51918, 74370,
                                                                       30028, 30128, 53928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 77268, 0, 3,
                                                                       74370, 52008, 74496,
                                                                       30128, 30228, 54078,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 77478, 0, 3,
                                                                       74496, 52098, 74622,
                                                                       30228, 30328, 54228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 77688, 0, 3,
                                                                       74748, 52278, 74958,
                                                                       30528, 30678, 54378,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 78003, 0, 3,
                                                                       74958, 52428, 75168,
                                                                       30678, 30828, 54603,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 78318, 0, 3,
                                                                       75168, 52578, 75378,
                                                                       30828, 30978, 54828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 78633, 0, 3,
                                                                       75378, 52728, 75588,
                                                                       30978, 31128, 55053,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 78948, 0, 3,
                                                                       75588, 52878, 75798,
                                                                       31128, 31278, 55278,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 79263, 0, 3,
                                                                       75798, 53028, 76008,
                                                                       31278, 31428, 55503,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 79578, 0, 3,
                                                                       76218, 53328, 76428,
                                                                       31728, 31878, 55728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 79893, 0, 3,
                                                                       76428, 53478, 76638,
                                                                       31878, 32028, 55953,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 80208, 0, 3,
                                                                       76638, 53628, 76848,
                                                                       32028, 32178, 56178,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 80523, 0, 3,
                                                                       76848, 53778, 77058,
                                                                       32178, 32328, 56403,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 80838, 0, 3,
                                                                       77058, 53928, 77268,
                                                                       32328, 32478, 56628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 81153, 0, 3,
                                                                       77268, 54078, 77478,
                                                                       32478, 32628, 56853,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 81468, 0, 3,
                                                                       77688, 54378, 78003,
                                                                       32928, 33138, 57078,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 81909, 0, 3,
                                                                       78003, 54603, 78318,
                                                                       33138, 33348, 57393,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 82350, 0, 3,
                                                                       78318, 54828, 78633,
                                                                       33348, 33558, 57708,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 82791, 0, 3,
                                                                       78633, 55053, 78948,
                                                                       33558, 33768, 58023,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 83232, 0, 3,
                                                                       78948, 55278, 79263,
                                                                       33768, 33978, 58338,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 83673, 0, 3,
                                                                       79578, 55728, 79893,
                                                                       34398, 34608, 58653,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 84114, 0, 3,
                                                                       79893, 55953, 80208,
                                                                       34608, 34818, 58968,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 84555, 0, 3,
                                                                       80208, 56178, 80523,
                                                                       34818, 35028, 59283,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 84996, 0, 3,
                                                                       80523, 56403, 80838,
                                                                       35028, 35238, 59598,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 85437, 0, 3,
                                                                       80838, 56628, 81153,
                                                                       35238, 35448, 59913,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 85878, 0, 3,
                                                                       81468, 57078, 81909,
                                                                       35868, 36148, 60228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 86466, 0, 3,
                                                                       81909, 57393, 82350,
                                                                       36148, 36428, 60648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 87054, 0, 3,
                                                                       82350, 57708, 82791,
                                                                       36428, 36708, 61068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 87642, 0, 3,
                                                                       82791, 58023, 83232,
                                                                       36708, 36988, 61488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 88230, 0, 3,
                                                                       83673, 58653, 84114,
                                                                       37548, 37828, 61908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 88818, 0, 3,
                                                                       84114, 58968, 84555,
                                                                       37828, 38108, 62328,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 89406, 0, 3,
                                                                       84555, 59283, 84996,
                                                                       38108, 38388, 62748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 89994, 0, 3,
                                                                       84996, 59598, 85437,
                                                                       38388, 38668, 63168,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 90582, 0, 3,
                                                                       85878, 60228, 86466,
                                                                       39228, 39588, 63588,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 91338, 0, 3,
                                                                       86466, 60648, 87054,
                                                                       39588, 39948, 64128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 92094, 0, 3,
                                                                       87054, 61068, 87642,
                                                                       39948, 40308, 64668,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 92850, 0, 3,
                                                                       88230, 61908, 88818,
                                                                       41028, 41388, 65208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 93606, 0, 3,
                                                                       88818, 62328, 89406,
                                                                       41388, 41748, 65748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 94362, 0, 3,
                                                                       89406, 62748, 89994,
                                                                       41748, 42108, 66288,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 95118, 0, 3,
                                                                       90582, 63588, 91338,
                                                                       42828, 43278, 66828,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 96063, 0, 3,
                                                                       91338, 64128, 92094,
                                                                       43278, 43728, 67503,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 97008, 0, 3,
                                                                       92850, 65208, 93606,
                                                                       44628, 45078, 68178,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 97953, 0, 3,
                                                                       93606, 65748, 94362,
                                                                       45078, 45528, 68853,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 98898, 0, 3,
                                                                       95118, 66828, 96063,
                                                                       46428, 46978, 69528,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 100053, 0, 3,
                                                                       97008, 68178, 97953,
                                                                       48078, 48628, 70353,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 101208, 81468, 441, ncols);

                    simdfunc::contract_primitives(buffer, 101880, 83673, 441, ncols);

                    simdfunc::contract_primitives(buffer, 102552, 85878, 588, ncols);

                    simdfunc::contract_primitives(buffer, 103448, 88230, 588, ncols);

                    simdfunc::contract_primitives(buffer, 104344, 90582, 756, ncols);

                    simdfunc::contract_primitives(buffer, 105496, 92850, 756, ncols);

                    simdfunc::contract_primitives(buffer, 106648, 95118, 945, ncols);

                    simdfunc::contract_primitives(buffer, 108088, 97008, 945, ncols);

                    simdfunc::contract_primitives(buffer, 109528, 98898, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 111288, 100053, 1155, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 101649, 101208, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 102321, 101880, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 103140, 102552, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 104036, 103448, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 105100, 104344, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 106252, 105496, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 107593, 106648, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 109033, 108088, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 110683, 109528, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 112443, 111288, 55, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 113048, 101649, 103140, 11,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 113741, 102321, 104036, 11,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 114434, 103140, 105100, 11,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 115358, 104036, 106252, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 116282, 105100, 107593, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 117470, 106252, 109033, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 118658, 107593, 110683, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 120143, 109033, 112443, 11,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 121628, 113048, 114434, 11,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 123014, 113741, 115358, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 124400, 114434, 116282, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 126248, 115358, 117470, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 128096, 116282, 118658, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 130472, 117470, 120143, 11,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 132848, 121628, 124400, 11,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 135158, 123014, 126248, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 137468, 124400, 128096, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 140548, 126248, 130472, 11,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 143628, 132848, 137468, 11,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 147093, 135158, 140548, 11,
                                             nmax);

        simdtrf::transform_g_inner(buffer, 150558, 147093, 21, 11, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 150558, 99, nmax);

        simdtrf::transform_g_inner(buffer, 150558, 143628, 21, 11, nmax);

        simdtrf::transform_h_outer(values + 1089 * nvalues + n * npairs, nvalues, buffer, 150558,
                                   99, nmax);
    }

    for (size_t m = 0; m < 2178; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
