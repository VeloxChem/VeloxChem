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


#include "SimdThreeCenterElectronRepulsionRsRecIGK.hpp"

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
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_igk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_igk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 425188, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3510 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 425188, 331768, 22470, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                                                            15, 16, 17}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 24, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16, 17}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 69, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 75, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 78, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 81, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 84, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 87, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 90, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 93, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 96, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 99, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 105, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 111, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 114, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 117, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 120, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 123, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 126, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 129, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 39, 40,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 135, 0, 3, 40, 41,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 7, 8,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 144, 0, 3, 8, 9,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 150, 0, 3, 9, 10,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 156, 0, 3, 10, 11,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 11, 12,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 12, 13,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 174, 0, 3, 13, 14,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 180, 0, 3, 14, 15,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 186, 0, 3, 15, 16,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 16, 17,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 17, 18,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 204, 0, 3, 18, 19,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 210, 0, 3, 19, 20,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 216, 0, 3, 20, 21,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 21, 22,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 25, 26,
                                                                       90, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 234, 0, 3, 26, 27,
                                                                       93, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 240, 0, 3, 27, 28,
                                                                       96, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 246, 0, 3, 28, 29,
                                                                       99, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 252, 0, 3, 29, 30,
                                                                       102, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 30, 31,
                                                                       105, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 264, 0, 3, 31, 32,
                                                                       108, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 270, 0, 3, 32, 33,
                                                                       111, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 276, 0, 3, 33, 34,
                                                                       114, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 282, 0, 3, 34, 35,
                                                                       117, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 35, 36,
                                                                       120, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 294, 0, 3, 36, 37,
                                                                       123, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 300, 0, 3, 37, 38,
                                                                       126, 129, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 306, 0, 3, 38, 39,
                                                                       129, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 312, 0, 3, 39, 40,
                                                                       132, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 42, 45,
                                                                       138, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 45, 48,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 48, 51,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 51, 54,
                                                                       156, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 54, 57,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 57, 60,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 60, 63,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 63, 66,
                                                                       180, 186, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 66, 69,
                                                                       186, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 69, 72,
                                                                       192, 198, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 72, 75,
                                                                       198, 204, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 75, 78,
                                                                       204, 210, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 78, 81,
                                                                       210, 216, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 81, 84,
                                                                       216, 222, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 90, 93,
                                                                       228, 234, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 93, 96,
                                                                       234, 240, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 96, 99,
                                                                       240, 246, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 99,
                                                                       102, 246, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 102,
                                                                       105, 252, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 105,
                                                                       108, 258, 264, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 108,
                                                                       111, 264, 270, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 528, 0, 3, 111,
                                                                       114, 270, 276, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 538, 0, 3, 114,
                                                                       117, 276, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 548, 0, 3, 117,
                                                                       120, 282, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 558, 0, 3, 120,
                                                                       123, 288, 294, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 568, 0, 3, 123,
                                                                       126, 294, 300, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 578, 0, 3, 126,
                                                                       129, 300, 306, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 588, 0, 3, 129,
                                                                       132, 306, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 598, 0, 3, 138,
                                                                       144, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 613, 0, 3, 144,
                                                                       150, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 628, 0, 3, 150,
                                                                       156, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 643, 0, 3, 156,
                                                                       162, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 658, 0, 3, 162,
                                                                       168, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 673, 0, 3, 168,
                                                                       174, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 688, 0, 3, 174,
                                                                       180, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 703, 0, 3, 180,
                                                                       186, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 718, 0, 3, 186,
                                                                       192, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 733, 0, 3, 192,
                                                                       198, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 748, 0, 3, 198,
                                                                       204, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 763, 0, 3, 204,
                                                                       210, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 778, 0, 3, 210,
                                                                       216, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 793, 0, 3, 228,
                                                                       234, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 808, 0, 3, 234,
                                                                       240, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 823, 0, 3, 240,
                                                                       246, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 838, 0, 3, 246,
                                                                       252, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 853, 0, 3, 252,
                                                                       258, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 868, 0, 3, 258,
                                                                       264, 508, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 883, 0, 3, 264,
                                                                       270, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 898, 0, 3, 270,
                                                                       276, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 913, 0, 3, 276,
                                                                       282, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 928, 0, 3, 282,
                                                                       288, 548, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 943, 0, 3, 288,
                                                                       294, 558, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 958, 0, 3, 294,
                                                                       300, 568, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 973, 0, 3, 300,
                                                                       306, 578, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 988, 0, 3, 318,
                                                                       328, 598, 613, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1009, 0, 3, 328,
                                                                       338, 613, 628, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 338,
                                                                       348, 628, 643, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1051, 0, 3, 348,
                                                                       358, 643, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 358,
                                                                       368, 658, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1093, 0, 3, 368,
                                                                       378, 673, 688, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 378,
                                                                       388, 688, 703, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1135, 0, 3, 388,
                                                                       398, 703, 718, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 398,
                                                                       408, 718, 733, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1177, 0, 3, 408,
                                                                       418, 733, 748, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1198, 0, 3, 418,
                                                                       428, 748, 763, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1219, 0, 3, 428,
                                                                       438, 763, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 458,
                                                                       468, 793, 808, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1261, 0, 3, 468,
                                                                       478, 808, 823, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1282, 0, 3, 478,
                                                                       488, 823, 838, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1303, 0, 3, 488,
                                                                       498, 838, 853, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 498,
                                                                       508, 853, 868, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1345, 0, 3, 508,
                                                                       518, 868, 883, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1366, 0, 3, 518,
                                                                       528, 883, 898, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1387, 0, 3, 528,
                                                                       538, 898, 913, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 538,
                                                                       548, 913, 928, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1429, 0, 3, 548,
                                                                       558, 928, 943, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1450, 0, 3, 558,
                                                                       568, 943, 958, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1471, 0, 3, 568,
                                                                       578, 958, 973, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 598,
                                                                       613, 988, 1009, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 613,
                                                                       628, 1009, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 628,
                                                                       643, 1030, 1051, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 643,
                                                                       658, 1051, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 658,
                                                                       673, 1072, 1093, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 673,
                                                                       688, 1093, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 688,
                                                                       703, 1114, 1135, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 703,
                                                                       718, 1135, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 718,
                                                                       733, 1156, 1177, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 733,
                                                                       748, 1177, 1198, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 748,
                                                                       763, 1198, 1219, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 793,
                                                                       808, 1240, 1261, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 808,
                                                                       823, 1261, 1282, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 823,
                                                                       838, 1282, 1303, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 838,
                                                                       853, 1303, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 853,
                                                                       868, 1324, 1345, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 868,
                                                                       883, 1345, 1366, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1968, 0, 3, 883,
                                                                       898, 1366, 1387, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1996, 0, 3, 898,
                                                                       913, 1387, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 913,
                                                                       928, 1408, 1429, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2052, 0, 3, 928,
                                                                       943, 1429, 1450, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 943,
                                                                       958, 1450, 1471, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 988,
                                                                       1009, 1492, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2144, 0, 3, 1009,
                                                                       1030, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2180, 0, 3, 1030,
                                                                       1051, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2216, 0, 3, 1051,
                                                                       1072, 1576, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2252, 0, 3, 1072,
                                                                       1093, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2288, 0, 3, 1093,
                                                                       1114, 1632, 1660, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2324, 0, 3, 1114,
                                                                       1135, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2360, 0, 3, 1135,
                                                                       1156, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2396, 0, 3, 1156,
                                                                       1177, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2432, 0, 3, 1177,
                                                                       1198, 1744, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2468, 0, 3, 1240,
                                                                       1261, 1800, 1828, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2504, 0, 3, 1261,
                                                                       1282, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2540, 0, 3, 1282,
                                                                       1303, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2576, 0, 3, 1303,
                                                                       1324, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2612, 0, 3, 1324,
                                                                       1345, 1912, 1940, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1345,
                                                                       1366, 1940, 1968, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2684, 0, 3, 1366,
                                                                       1387, 1968, 1996, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2720, 0, 3, 1387,
                                                                       1408, 1996, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2756, 0, 3, 1408,
                                                                       1429, 2024, 2052, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2792, 0, 3, 1429,
                                                                       1450, 2052, 2080, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2828, 0, 3, 1492,
                                                                       1520, 2108, 2144, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2873, 0, 3, 1520,
                                                                       1548, 2144, 2180, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2918, 0, 3, 1548,
                                                                       1576, 2180, 2216, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2963, 0, 3, 1576,
                                                                       1604, 2216, 2252, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3008, 0, 3, 1604,
                                                                       1632, 2252, 2288, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3053, 0, 3, 1632,
                                                                       1660, 2288, 2324, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3098, 0, 3, 1660,
                                                                       1688, 2324, 2360, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3143, 0, 3, 1688,
                                                                       1716, 2360, 2396, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3188, 0, 3, 1716,
                                                                       1744, 2396, 2432, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3233, 0, 3, 1800,
                                                                       1828, 2468, 2504, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3278, 0, 3, 1828,
                                                                       1856, 2504, 2540, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3323, 0, 3, 1856,
                                                                       1884, 2540, 2576, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3368, 0, 3, 1884,
                                                                       1912, 2576, 2612, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3413, 0, 3, 1912,
                                                                       1940, 2612, 2648, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3458, 0, 3, 1940,
                                                                       1968, 2648, 2684, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3503, 0, 3, 1968,
                                                                       1996, 2684, 2720, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3548, 0, 3, 1996,
                                                                       2024, 2720, 2756, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3593, 0, 3, 2024,
                                                                       2052, 2756, 2792, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3638, 0, 3, 2108,
                                                                       2144, 2828, 2873, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3693, 0, 3, 2144,
                                                                       2180, 2873, 2918, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3748, 0, 3, 2180,
                                                                       2216, 2918, 2963, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3803, 0, 3, 2216,
                                                                       2252, 2963, 3008, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3858, 0, 3, 2252,
                                                                       2288, 3008, 3053, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3913, 0, 3, 2288,
                                                                       2324, 3053, 3098, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 2324,
                                                                       2360, 3098, 3143, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4023, 0, 3, 2360,
                                                                       2396, 3143, 3188, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4078, 0, 3, 2468,
                                                                       2504, 3233, 3278, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4133, 0, 3, 2504,
                                                                       2540, 3278, 3323, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4188, 0, 3, 2540,
                                                                       2576, 3323, 3368, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4243, 0, 3, 2576,
                                                                       2612, 3368, 3413, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 2612,
                                                                       2648, 3413, 3458, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4353, 0, 3, 2648,
                                                                       2684, 3458, 3503, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4408, 0, 3, 2684,
                                                                       2720, 3503, 3548, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4463, 0, 3, 2720,
                                                                       2756, 3548, 3593, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4518, 0, 3, 2828,
                                                                       2873, 3638, 3693, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4584, 0, 3, 2873,
                                                                       2918, 3693, 3748, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4650, 0, 3, 2918,
                                                                       2963, 3748, 3803, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4716, 0, 3, 2963,
                                                                       3008, 3803, 3858, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4782, 0, 3, 3008,
                                                                       3053, 3858, 3913, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4848, 0, 3, 3053,
                                                                       3098, 3913, 3968, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4914, 0, 3, 3098,
                                                                       3143, 3968, 4023, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4980, 0, 3, 3233,
                                                                       3278, 4078, 4133, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5046, 0, 3, 3278,
                                                                       3323, 4133, 4188, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5112, 0, 3, 3323,
                                                                       3368, 4188, 4243, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5178, 0, 3, 3368,
                                                                       3413, 4243, 4298, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5244, 0, 3, 3413,
                                                                       3458, 4298, 4353, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5310, 0, 3, 3458,
                                                                       3503, 4353, 4408, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5376, 0, 3, 3503,
                                                                       3548, 4408, 4463, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5442, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5445, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5448, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5451, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5454, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5457, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5460, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5463, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5466, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5469, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5472, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5475, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5478, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5481, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5484, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5487, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5490, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5493, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5496, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5499, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5502, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5505, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5508, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5511, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5514, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5517, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5520, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5523, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5526, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5529, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5532, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5535, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5538, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5541, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5544, 3, 7, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5553, 3, 8, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5562, 3, 9, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5571, 3, 10, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5580, 3, 11, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5589, 3, 12, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5598, 3, 13, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5607, 3, 14, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5616, 3, 15, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5625, 3, 16, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5634, 3, 17, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5643, 3, 18, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5652, 3, 19, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5661, 3, 20, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5670, 3, 21, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5679, 3, 22, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5688, 3, 25, 90,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5697, 3, 26, 93,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5706, 3, 27, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5715, 3, 28, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5724, 3, 29, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5733, 3, 30, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5742, 3, 31, 108,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5751, 3, 32, 111,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5760, 3, 33, 114,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5769, 3, 34, 117,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5778, 3, 35, 120,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5787, 3, 36, 123,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5796, 3, 37, 126,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5805, 3, 38, 129,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5814, 3, 39, 132,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5823, 3, 40, 135,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5832, 3, 42, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5850, 3, 45, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5868, 3, 48, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5886, 3, 51, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5904, 3, 54, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5922, 3, 57, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5940, 3, 60, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5958, 3, 63, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5976, 3, 66, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5994, 3, 69, 192,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6012, 3, 72, 198,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6030, 3, 75, 204,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6048, 3, 78, 210,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6066, 3, 81, 216,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6084, 3, 84, 222,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6102, 3, 90, 228,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6120, 3, 93, 234,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6138, 3, 96, 240,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6156, 3, 99, 246,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6174, 3, 102, 252,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6192, 3, 105, 258,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6210, 3, 108, 264,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6228, 3, 111, 270,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6246, 3, 114, 276,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6264, 3, 117, 282,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6282, 3, 120, 288,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6300, 3, 123, 294,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6318, 3, 126, 300,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6336, 3, 129, 306,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6354, 3, 132, 312,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6372, 3, 138, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6402, 3, 144, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6432, 3, 150, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6462, 3, 156, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6492, 3, 162, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6522, 3, 168, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6552, 3, 174, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6582, 3, 180, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6612, 3, 186, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6642, 3, 192, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6672, 3, 198, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6702, 3, 204, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6732, 3, 210, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6762, 3, 216, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6792, 3, 228, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6822, 3, 234, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6852, 3, 240, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6882, 3, 246, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6912, 3, 252, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6942, 3, 258, 508,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6972, 3, 264, 518,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7002, 3, 270, 528,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7032, 3, 276, 538,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7062, 3, 282, 548,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7092, 3, 288, 558,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7122, 3, 294, 568,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7152, 3, 300, 578,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7182, 3, 306, 588,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7212, 3, 318, 598,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7257, 3, 328, 613,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7302, 3, 338, 628,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7347, 3, 348, 643,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7392, 3, 358, 658,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7437, 3, 368, 673,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7482, 3, 378, 688,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7527, 3, 388, 703,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7572, 3, 398, 718,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7617, 3, 408, 733,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7662, 3, 418, 748,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7707, 3, 428, 763,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7752, 3, 438, 778,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7797, 3, 458, 793,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7842, 3, 468, 808,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7887, 3, 478, 823,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7932, 3, 488, 838,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7977, 3, 498, 853,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8022, 3, 508, 868,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8067, 3, 518, 883,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8112, 3, 528, 898,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8157, 3, 538, 913,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8202, 3, 548, 928,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8247, 3, 558, 943,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8292, 3, 568, 958,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8337, 3, 578, 973,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8382, 3, 598, 988,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8445, 3, 613,
                                                                       1009, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8508, 3, 628,
                                                                       1030, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8571, 3, 643,
                                                                       1051, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8634, 3, 658,
                                                                       1072, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8697, 3, 673,
                                                                       1093, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8760, 3, 688,
                                                                       1114, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8823, 3, 703,
                                                                       1135, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8886, 3, 718,
                                                                       1156, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8949, 3, 733,
                                                                       1177, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9012, 3, 748,
                                                                       1198, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9075, 3, 763,
                                                                       1219, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9138, 3, 793,
                                                                       1240, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9201, 3, 808,
                                                                       1261, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9264, 3, 823,
                                                                       1282, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9327, 3, 838,
                                                                       1303, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9390, 3, 853,
                                                                       1324, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9453, 3, 868,
                                                                       1345, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9516, 3, 883,
                                                                       1366, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9579, 3, 898,
                                                                       1387, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9642, 3, 913,
                                                                       1408, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9705, 3, 928,
                                                                       1429, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9768, 3, 943,
                                                                       1450, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9831, 3, 958,
                                                                       1471, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9894, 3, 988,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9978, 3, 1009,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10062, 3, 1030,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10146, 3, 1051,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10230, 3, 1072,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10314, 3, 1093,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10398, 3, 1114,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10482, 3, 1135,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10566, 3, 1156,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10650, 3, 1177,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10734, 3, 1198,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10818, 3, 1240,
                                                                       1800, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10902, 3, 1261,
                                                                       1828, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10986, 3, 1282,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11070, 3, 1303,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11154, 3, 1324,
                                                                       1912, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11238, 3, 1345,
                                                                       1940, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11322, 3, 1366,
                                                                       1968, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11406, 3, 1387,
                                                                       1996, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11490, 3, 1408,
                                                                       2024, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11574, 3, 1429,
                                                                       2052, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11658, 3, 1450,
                                                                       2080, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11742, 3, 1492,
                                                                       2108, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11850, 3, 1520,
                                                                       2144, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11958, 3, 1548,
                                                                       2180, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12066, 3, 1576,
                                                                       2216, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12174, 3, 1604,
                                                                       2252, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12282, 3, 1632,
                                                                       2288, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12390, 3, 1660,
                                                                       2324, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12498, 3, 1688,
                                                                       2360, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12606, 3, 1716,
                                                                       2396, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12714, 3, 1744,
                                                                       2432, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12822, 3, 1800,
                                                                       2468, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12930, 3, 1828,
                                                                       2504, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13038, 3, 1856,
                                                                       2540, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13146, 3, 1884,
                                                                       2576, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13254, 3, 1912,
                                                                       2612, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13362, 3, 1940,
                                                                       2648, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13470, 3, 1968,
                                                                       2684, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13578, 3, 1996,
                                                                       2720, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13686, 3, 2024,
                                                                       2756, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13794, 3, 2052,
                                                                       2792, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13902, 3, 2108,
                                                                       2828, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14037, 3, 2144,
                                                                       2873, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14172, 3, 2180,
                                                                       2918, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14307, 3, 2216,
                                                                       2963, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14442, 3, 2252,
                                                                       3008, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14577, 3, 2288,
                                                                       3053, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14712, 3, 2324,
                                                                       3098, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14847, 3, 2360,
                                                                       3143, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14982, 3, 2396,
                                                                       3188, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15117, 3, 2468,
                                                                       3233, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15252, 3, 2504,
                                                                       3278, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15387, 3, 2540,
                                                                       3323, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15522, 3, 2576,
                                                                       3368, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15657, 3, 2612,
                                                                       3413, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15792, 3, 2648,
                                                                       3458, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15927, 3, 2684,
                                                                       3503, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16062, 3, 2720,
                                                                       3548, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16197, 3, 2756,
                                                                       3593, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16332, 3, 2828,
                                                                       3638, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16497, 3, 2873,
                                                                       3693, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16662, 3, 2918,
                                                                       3748, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16827, 3, 2963,
                                                                       3803, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16992, 3, 3008,
                                                                       3858, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17157, 3, 3053,
                                                                       3913, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17322, 3, 3098,
                                                                       3968, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17487, 3, 3143,
                                                                       4023, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17652, 3, 3233,
                                                                       4078, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17817, 3, 3278,
                                                                       4133, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17982, 3, 3323,
                                                                       4188, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18147, 3, 3368,
                                                                       4243, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18312, 3, 3413,
                                                                       4298, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18477, 3, 3458,
                                                                       4353, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18642, 3, 3503,
                                                                       4408, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18807, 3, 3548,
                                                                       4463, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18972, 3, 3638,
                                                                       4518, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19170, 3, 3693,
                                                                       4584, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19368, 3, 3748,
                                                                       4650, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19566, 3, 3803,
                                                                       4716, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19764, 3, 3858,
                                                                       4782, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19962, 3, 3913,
                                                                       4848, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20160, 3, 3968,
                                                                       4914, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20358, 3, 4078,
                                                                       4980, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20556, 3, 4133,
                                                                       5046, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20754, 3, 4188,
                                                                       5112, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20952, 3, 4243,
                                                                       5178, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21150, 3, 4298,
                                                                       5244, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21348, 3, 4353,
                                                                       5310, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21546, 3, 4408,
                                                                       5376, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21744, 3, 7, 8,
                                                                       5448, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21750, 3, 8, 9,
                                                                       5451, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21756, 3, 9, 10,
                                                                       5454, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21762, 3, 10, 11,
                                                                       5457, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21768, 3, 11, 12,
                                                                       5460, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21774, 3, 12, 13,
                                                                       5463, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21780, 3, 13, 14,
                                                                       5466, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21786, 3, 14, 15,
                                                                       5469, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21792, 3, 15, 16,
                                                                       5472, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21798, 3, 16, 17,
                                                                       5475, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21804, 3, 17, 18,
                                                                       5478, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21810, 3, 18, 19,
                                                                       5481, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21816, 3, 19, 20,
                                                                       5484, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21822, 3, 20, 21,
                                                                       5487, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21828, 3, 21, 22,
                                                                       5490, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21834, 3, 25, 26,
                                                                       5499, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21840, 3, 26, 27,
                                                                       5502, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21846, 3, 27, 28,
                                                                       5505, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21852, 3, 28, 29,
                                                                       5508, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21858, 3, 29, 30,
                                                                       5511, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21864, 3, 30, 31,
                                                                       5514, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21870, 3, 31, 32,
                                                                       5517, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21876, 3, 32, 33,
                                                                       5520, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21882, 3, 33, 34,
                                                                       5523, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21888, 3, 34, 35,
                                                                       5526, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21894, 3, 35, 36,
                                                                       5529, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21900, 3, 36, 37,
                                                                       5532, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21906, 3, 37, 38,
                                                                       5535, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21912, 3, 38, 39,
                                                                       5538, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21918, 3, 39, 40,
                                                                       5541, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 21924, 0, 3,
                                                                       21744, 5448, 21750, 5562,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 21942, 0, 3,
                                                                       21750, 5451, 21756, 5571,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 21960, 0, 3,
                                                                       21756, 5454, 21762, 5580,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 21978, 0, 3,
                                                                       21762, 5457, 21768, 5589,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 21996, 0, 3,
                                                                       21768, 5460, 21774, 5598,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22014, 0, 3,
                                                                       21774, 5463, 21780, 5607,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22032, 0, 3,
                                                                       21780, 5466, 21786, 5616,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22050, 0, 3,
                                                                       21786, 5469, 21792, 5625,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22068, 0, 3,
                                                                       21792, 5472, 21798, 5634,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22086, 0, 3,
                                                                       21798, 5475, 21804, 5643,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22104, 0, 3,
                                                                       21804, 5478, 21810, 5652,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22122, 0, 3,
                                                                       21810, 5481, 21816, 5661,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22140, 0, 3,
                                                                       21816, 5484, 21822, 5670,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22158, 0, 3,
                                                                       21822, 5487, 21828, 5679,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22176, 0, 3,
                                                                       21834, 5499, 21840, 5706,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22194, 0, 3,
                                                                       21840, 5502, 21846, 5715,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22212, 0, 3,
                                                                       21846, 5505, 21852, 5724,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22230, 0, 3,
                                                                       21852, 5508, 21858, 5733,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22248, 0, 3,
                                                                       21858, 5511, 21864, 5742,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22266, 0, 3,
                                                                       21864, 5514, 21870, 5751,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22284, 0, 3,
                                                                       21870, 5517, 21876, 5760,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22302, 0, 3,
                                                                       21876, 5520, 21882, 5769,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22320, 0, 3,
                                                                       21882, 5523, 21888, 5778,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22338, 0, 3,
                                                                       21888, 5526, 21894, 5787,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22356, 0, 3,
                                                                       21894, 5529, 21900, 5796,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22374, 0, 3,
                                                                       21900, 5532, 21906, 5805,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22392, 0, 3,
                                                                       21906, 5535, 21912, 5814,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22410, 0, 3,
                                                                       21912, 5538, 21918, 5823,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22428, 0, 3,
                                                                       21924, 5562, 21942, 138,
                                                                       144, 5868, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22464, 0, 3,
                                                                       21942, 5571, 21960, 144,
                                                                       150, 5886, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22500, 0, 3,
                                                                       21960, 5580, 21978, 150,
                                                                       156, 5904, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22536, 0, 3,
                                                                       21978, 5589, 21996, 156,
                                                                       162, 5922, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22572, 0, 3,
                                                                       21996, 5598, 22014, 162,
                                                                       168, 5940, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22608, 0, 3,
                                                                       22014, 5607, 22032, 168,
                                                                       174, 5958, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22644, 0, 3,
                                                                       22032, 5616, 22050, 174,
                                                                       180, 5976, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22680, 0, 3,
                                                                       22050, 5625, 22068, 180,
                                                                       186, 5994, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22716, 0, 3,
                                                                       22068, 5634, 22086, 186,
                                                                       192, 6012, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22752, 0, 3,
                                                                       22086, 5643, 22104, 192,
                                                                       198, 6030, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22788, 0, 3,
                                                                       22104, 5652, 22122, 198,
                                                                       204, 6048, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22824, 0, 3,
                                                                       22122, 5661, 22140, 204,
                                                                       210, 6066, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22860, 0, 3,
                                                                       22140, 5670, 22158, 210,
                                                                       216, 6084, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22896, 0, 3,
                                                                       22176, 5706, 22194, 228,
                                                                       234, 6138, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22932, 0, 3,
                                                                       22194, 5715, 22212, 234,
                                                                       240, 6156, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22968, 0, 3,
                                                                       22212, 5724, 22230, 240,
                                                                       246, 6174, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23004, 0, 3,
                                                                       22230, 5733, 22248, 246,
                                                                       252, 6192, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23040, 0, 3,
                                                                       22248, 5742, 22266, 252,
                                                                       258, 6210, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23076, 0, 3,
                                                                       22266, 5751, 22284, 258,
                                                                       264, 6228, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23112, 0, 3,
                                                                       22284, 5760, 22302, 264,
                                                                       270, 6246, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23148, 0, 3,
                                                                       22302, 5769, 22320, 270,
                                                                       276, 6264, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23184, 0, 3,
                                                                       22320, 5778, 22338, 276,
                                                                       282, 6282, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23220, 0, 3,
                                                                       22338, 5787, 22356, 282,
                                                                       288, 6300, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23256, 0, 3,
                                                                       22356, 5796, 22374, 288,
                                                                       294, 6318, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23292, 0, 3,
                                                                       22374, 5805, 22392, 294,
                                                                       300, 6336, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23328, 0, 3,
                                                                       22392, 5814, 22410, 300,
                                                                       306, 6354, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23364, 0, 3,
                                                                       22428, 5868, 22464, 318,
                                                                       328, 6432, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23424, 0, 3,
                                                                       22464, 5886, 22500, 328,
                                                                       338, 6462, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23484, 0, 3,
                                                                       22500, 5904, 22536, 338,
                                                                       348, 6492, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23544, 0, 3,
                                                                       22536, 5922, 22572, 348,
                                                                       358, 6522, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23604, 0, 3,
                                                                       22572, 5940, 22608, 358,
                                                                       368, 6552, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23664, 0, 3,
                                                                       22608, 5958, 22644, 368,
                                                                       378, 6582, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23724, 0, 3,
                                                                       22644, 5976, 22680, 378,
                                                                       388, 6612, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23784, 0, 3,
                                                                       22680, 5994, 22716, 388,
                                                                       398, 6642, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23844, 0, 3,
                                                                       22716, 6012, 22752, 398,
                                                                       408, 6672, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23904, 0, 3,
                                                                       22752, 6030, 22788, 408,
                                                                       418, 6702, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23964, 0, 3,
                                                                       22788, 6048, 22824, 418,
                                                                       428, 6732, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24024, 0, 3,
                                                                       22824, 6066, 22860, 428,
                                                                       438, 6762, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24084, 0, 3,
                                                                       22896, 6138, 22932, 458,
                                                                       468, 6852, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24144, 0, 3,
                                                                       22932, 6156, 22968, 468,
                                                                       478, 6882, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24204, 0, 3,
                                                                       22968, 6174, 23004, 478,
                                                                       488, 6912, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24264, 0, 3,
                                                                       23004, 6192, 23040, 488,
                                                                       498, 6942, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24324, 0, 3,
                                                                       23040, 6210, 23076, 498,
                                                                       508, 6972, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24384, 0, 3,
                                                                       23076, 6228, 23112, 508,
                                                                       518, 7002, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24444, 0, 3,
                                                                       23112, 6246, 23148, 518,
                                                                       528, 7032, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24504, 0, 3,
                                                                       23148, 6264, 23184, 528,
                                                                       538, 7062, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24564, 0, 3,
                                                                       23184, 6282, 23220, 538,
                                                                       548, 7092, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24624, 0, 3,
                                                                       23220, 6300, 23256, 548,
                                                                       558, 7122, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24684, 0, 3,
                                                                       23256, 6318, 23292, 558,
                                                                       568, 7152, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24744, 0, 3,
                                                                       23292, 6336, 23328, 568,
                                                                       578, 7182, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 24804, 0, 3,
                                                                       23364, 6432, 23424, 598,
                                                                       613, 7302, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 24894, 0, 3,
                                                                       23424, 6462, 23484, 613,
                                                                       628, 7347, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 24984, 0, 3,
                                                                       23484, 6492, 23544, 628,
                                                                       643, 7392, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25074, 0, 3,
                                                                       23544, 6522, 23604, 643,
                                                                       658, 7437, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25164, 0, 3,
                                                                       23604, 6552, 23664, 658,
                                                                       673, 7482, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25254, 0, 3,
                                                                       23664, 6582, 23724, 673,
                                                                       688, 7527, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25344, 0, 3,
                                                                       23724, 6612, 23784, 688,
                                                                       703, 7572, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25434, 0, 3,
                                                                       23784, 6642, 23844, 703,
                                                                       718, 7617, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25524, 0, 3,
                                                                       23844, 6672, 23904, 718,
                                                                       733, 7662, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25614, 0, 3,
                                                                       23904, 6702, 23964, 733,
                                                                       748, 7707, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25704, 0, 3,
                                                                       23964, 6732, 24024, 748,
                                                                       763, 7752, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25794, 0, 3,
                                                                       24084, 6852, 24144, 793,
                                                                       808, 7887, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25884, 0, 3,
                                                                       24144, 6882, 24204, 808,
                                                                       823, 7932, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25974, 0, 3,
                                                                       24204, 6912, 24264, 823,
                                                                       838, 7977, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26064, 0, 3,
                                                                       24264, 6942, 24324, 838,
                                                                       853, 8022, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26154, 0, 3,
                                                                       24324, 6972, 24384, 853,
                                                                       868, 8067, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26244, 0, 3,
                                                                       24384, 7002, 24444, 868,
                                                                       883, 8112, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26334, 0, 3,
                                                                       24444, 7032, 24504, 883,
                                                                       898, 8157, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26424, 0, 3,
                                                                       24504, 7062, 24564, 898,
                                                                       913, 8202, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26514, 0, 3,
                                                                       24564, 7092, 24624, 913,
                                                                       928, 8247, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26604, 0, 3,
                                                                       24624, 7122, 24684, 928,
                                                                       943, 8292, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26694, 0, 3,
                                                                       24684, 7152, 24744, 943,
                                                                       958, 8337, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 26784, 0, 3,
                                                                       24804, 7302, 24894, 988,
                                                                       1009, 8508, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 26910, 0, 3,
                                                                       24894, 7347, 24984, 1009,
                                                                       1030, 8571, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27036, 0, 3,
                                                                       24984, 7392, 25074, 1030,
                                                                       1051, 8634, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27162, 0, 3,
                                                                       25074, 7437, 25164, 1051,
                                                                       1072, 8697, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27288, 0, 3,
                                                                       25164, 7482, 25254, 1072,
                                                                       1093, 8760, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27414, 0, 3,
                                                                       25254, 7527, 25344, 1093,
                                                                       1114, 8823, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27540, 0, 3,
                                                                       25344, 7572, 25434, 1114,
                                                                       1135, 8886, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27666, 0, 3,
                                                                       25434, 7617, 25524, 1135,
                                                                       1156, 8949, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27792, 0, 3,
                                                                       25524, 7662, 25614, 1156,
                                                                       1177, 9012, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27918, 0, 3,
                                                                       25614, 7707, 25704, 1177,
                                                                       1198, 9075, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28044, 0, 3,
                                                                       25794, 7887, 25884, 1240,
                                                                       1261, 9264, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28170, 0, 3,
                                                                       25884, 7932, 25974, 1261,
                                                                       1282, 9327, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28296, 0, 3,
                                                                       25974, 7977, 26064, 1282,
                                                                       1303, 9390, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28422, 0, 3,
                                                                       26064, 8022, 26154, 1303,
                                                                       1324, 9453, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28548, 0, 3,
                                                                       26154, 8067, 26244, 1324,
                                                                       1345, 9516, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28674, 0, 3,
                                                                       26244, 8112, 26334, 1345,
                                                                       1366, 9579, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28800, 0, 3,
                                                                       26334, 8157, 26424, 1366,
                                                                       1387, 9642, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28926, 0, 3,
                                                                       26424, 8202, 26514, 1387,
                                                                       1408, 9705, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29052, 0, 3,
                                                                       26514, 8247, 26604, 1408,
                                                                       1429, 9768, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29178, 0, 3,
                                                                       26604, 8292, 26694, 1429,
                                                                       1450, 9831, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29304, 0, 3,
                                                                       26784, 8508, 26910, 1492,
                                                                       1520, 10062, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29472, 0, 3,
                                                                       26910, 8571, 27036, 1520,
                                                                       1548, 10146, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29640, 0, 3,
                                                                       27036, 8634, 27162, 1548,
                                                                       1576, 10230, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29808, 0, 3,
                                                                       27162, 8697, 27288, 1576,
                                                                       1604, 10314, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29976, 0, 3,
                                                                       27288, 8760, 27414, 1604,
                                                                       1632, 10398, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30144, 0, 3,
                                                                       27414, 8823, 27540, 1632,
                                                                       1660, 10482, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30312, 0, 3,
                                                                       27540, 8886, 27666, 1660,
                                                                       1688, 10566, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30480, 0, 3,
                                                                       27666, 8949, 27792, 1688,
                                                                       1716, 10650, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30648, 0, 3,
                                                                       27792, 9012, 27918, 1716,
                                                                       1744, 10734, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30816, 0, 3,
                                                                       28044, 9264, 28170, 1800,
                                                                       1828, 10986, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30984, 0, 3,
                                                                       28170, 9327, 28296, 1828,
                                                                       1856, 11070, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31152, 0, 3,
                                                                       28296, 9390, 28422, 1856,
                                                                       1884, 11154, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31320, 0, 3,
                                                                       28422, 9453, 28548, 1884,
                                                                       1912, 11238, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31488, 0, 3,
                                                                       28548, 9516, 28674, 1912,
                                                                       1940, 11322, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31656, 0, 3,
                                                                       28674, 9579, 28800, 1940,
                                                                       1968, 11406, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31824, 0, 3,
                                                                       28800, 9642, 28926, 1968,
                                                                       1996, 11490, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31992, 0, 3,
                                                                       28926, 9705, 29052, 1996,
                                                                       2024, 11574, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32160, 0, 3,
                                                                       29052, 9768, 29178, 2024,
                                                                       2052, 11658, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32328, 0, 3,
                                                                       29304, 10062, 29472, 2108,
                                                                       2144, 11958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32544, 0, 3,
                                                                       29472, 10146, 29640, 2144,
                                                                       2180, 12066, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32760, 0, 3,
                                                                       29640, 10230, 29808, 2180,
                                                                       2216, 12174, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32976, 0, 3,
                                                                       29808, 10314, 29976, 2216,
                                                                       2252, 12282, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 33192, 0, 3,
                                                                       29976, 10398, 30144, 2252,
                                                                       2288, 12390, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 33408, 0, 3,
                                                                       30144, 10482, 30312, 2288,
                                                                       2324, 12498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 33624, 0, 3,
                                                                       30312, 10566, 30480, 2324,
                                                                       2360, 12606, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 33840, 0, 3,
                                                                       30480, 10650, 30648, 2360,
                                                                       2396, 12714, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34056, 0, 3,
                                                                       30816, 10986, 30984, 2468,
                                                                       2504, 13038, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34272, 0, 3,
                                                                       30984, 11070, 31152, 2504,
                                                                       2540, 13146, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34488, 0, 3,
                                                                       31152, 11154, 31320, 2540,
                                                                       2576, 13254, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34704, 0, 3,
                                                                       31320, 11238, 31488, 2576,
                                                                       2612, 13362, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34920, 0, 3,
                                                                       31488, 11322, 31656, 2612,
                                                                       2648, 13470, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 35136, 0, 3,
                                                                       31656, 11406, 31824, 2648,
                                                                       2684, 13578, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 35352, 0, 3,
                                                                       31824, 11490, 31992, 2684,
                                                                       2720, 13686, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 35568, 0, 3,
                                                                       31992, 11574, 32160, 2720,
                                                                       2756, 13794, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35784, 0, 3,
                                                                       32328, 11958, 32544, 2828,
                                                                       2873, 14172, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 36054, 0, 3,
                                                                       32544, 12066, 32760, 2873,
                                                                       2918, 14307, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 36324, 0, 3,
                                                                       32760, 12174, 32976, 2918,
                                                                       2963, 14442, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 36594, 0, 3,
                                                                       32976, 12282, 33192, 2963,
                                                                       3008, 14577, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 36864, 0, 3,
                                                                       33192, 12390, 33408, 3008,
                                                                       3053, 14712, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 37134, 0, 3,
                                                                       33408, 12498, 33624, 3053,
                                                                       3098, 14847, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 37404, 0, 3,
                                                                       33624, 12606, 33840, 3098,
                                                                       3143, 14982, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 37674, 0, 3,
                                                                       34056, 13038, 34272, 3233,
                                                                       3278, 15387, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 37944, 0, 3,
                                                                       34272, 13146, 34488, 3278,
                                                                       3323, 15522, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 38214, 0, 3,
                                                                       34488, 13254, 34704, 3323,
                                                                       3368, 15657, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 38484, 0, 3,
                                                                       34704, 13362, 34920, 3368,
                                                                       3413, 15792, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 38754, 0, 3,
                                                                       34920, 13470, 35136, 3413,
                                                                       3458, 15927, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 39024, 0, 3,
                                                                       35136, 13578, 35352, 3458,
                                                                       3503, 16062, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 39294, 0, 3,
                                                                       35352, 13686, 35568, 3503,
                                                                       3548, 16197, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 39564, 0, 3,
                                                                       35784, 14172, 36054, 3638,
                                                                       3693, 16662, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 39894, 0, 3,
                                                                       36054, 14307, 36324, 3693,
                                                                       3748, 16827, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 40224, 0, 3,
                                                                       36324, 14442, 36594, 3748,
                                                                       3803, 16992, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 40554, 0, 3,
                                                                       36594, 14577, 36864, 3803,
                                                                       3858, 17157, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 40884, 0, 3,
                                                                       36864, 14712, 37134, 3858,
                                                                       3913, 17322, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 41214, 0, 3,
                                                                       37134, 14847, 37404, 3913,
                                                                       3968, 17487, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 41544, 0, 3,
                                                                       37674, 15387, 37944, 4078,
                                                                       4133, 17982, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 41874, 0, 3,
                                                                       37944, 15522, 38214, 4133,
                                                                       4188, 18147, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 42204, 0, 3,
                                                                       38214, 15657, 38484, 4188,
                                                                       4243, 18312, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 42534, 0, 3,
                                                                       38484, 15792, 38754, 4243,
                                                                       4298, 18477, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 42864, 0, 3,
                                                                       38754, 15927, 39024, 4298,
                                                                       4353, 18642, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 43194, 0, 3,
                                                                       39024, 16062, 39294, 4353,
                                                                       4408, 18807, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 43524, 0, 3,
                                                                       39564, 16662, 39894, 4518,
                                                                       4584, 19368, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 43920, 0, 3,
                                                                       39894, 16827, 40224, 4584,
                                                                       4650, 19566, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 44316, 0, 3,
                                                                       40224, 16992, 40554, 4650,
                                                                       4716, 19764, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 44712, 0, 3,
                                                                       40554, 17157, 40884, 4716,
                                                                       4782, 19962, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 45108, 0, 3,
                                                                       40884, 17322, 41214, 4782,
                                                                       4848, 20160, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 45504, 0, 3,
                                                                       41544, 17982, 41874, 4980,
                                                                       5046, 20754, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 45900, 0, 3,
                                                                       41874, 18147, 42204, 5046,
                                                                       5112, 20952, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 46296, 0, 3,
                                                                       42204, 18312, 42534, 5112,
                                                                       5178, 21150, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 46692, 0, 3,
                                                                       42534, 18477, 42864, 5178,
                                                                       5244, 21348, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 47088, 0, 3,
                                                                       42864, 18642, 43194, 5244,
                                                                       5310, 21546, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47484, 3, 5442,
                                                                       5445, 21744, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47494, 3, 5445,
                                                                       5448, 21750, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47504, 3, 5448,
                                                                       5451, 21756, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47514, 3, 5451,
                                                                       5454, 21762, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47524, 3, 5454,
                                                                       5457, 21768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47534, 3, 5457,
                                                                       5460, 21774, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47544, 3, 5460,
                                                                       5463, 21780, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47554, 3, 5463,
                                                                       5466, 21786, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47564, 3, 5466,
                                                                       5469, 21792, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47574, 3, 5469,
                                                                       5472, 21798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47584, 3, 5472,
                                                                       5475, 21804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47594, 3, 5475,
                                                                       5478, 21810, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47604, 3, 5478,
                                                                       5481, 21816, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47614, 3, 5481,
                                                                       5484, 21822, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47624, 3, 5484,
                                                                       5487, 21828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47634, 3, 5493,
                                                                       5496, 21834, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47644, 3, 5496,
                                                                       5499, 21840, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47654, 3, 5499,
                                                                       5502, 21846, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47664, 3, 5502,
                                                                       5505, 21852, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47674, 3, 5505,
                                                                       5508, 21858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47684, 3, 5508,
                                                                       5511, 21864, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47694, 3, 5511,
                                                                       5514, 21870, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47704, 3, 5514,
                                                                       5517, 21876, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47714, 3, 5517,
                                                                       5520, 21882, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47724, 3, 5520,
                                                                       5523, 21888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47734, 3, 5523,
                                                                       5526, 21894, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47744, 3, 5526,
                                                                       5529, 21900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47754, 3, 5529,
                                                                       5532, 21906, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47764, 3, 5532,
                                                                       5535, 21912, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47774, 3, 5535,
                                                                       5538, 21918, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47784, 0, 3,
                                                                       47484, 21744, 47494, 5544,
                                                                       5553, 21924, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47814, 0, 3,
                                                                       47494, 21750, 47504, 5553,
                                                                       5562, 21942, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47844, 0, 3,
                                                                       47504, 21756, 47514, 5562,
                                                                       5571, 21960, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47874, 0, 3,
                                                                       47514, 21762, 47524, 5571,
                                                                       5580, 21978, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47904, 0, 3,
                                                                       47524, 21768, 47534, 5580,
                                                                       5589, 21996, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47934, 0, 3,
                                                                       47534, 21774, 47544, 5589,
                                                                       5598, 22014, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47964, 0, 3,
                                                                       47544, 21780, 47554, 5598,
                                                                       5607, 22032, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47994, 0, 3,
                                                                       47554, 21786, 47564, 5607,
                                                                       5616, 22050, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48024, 0, 3,
                                                                       47564, 21792, 47574, 5616,
                                                                       5625, 22068, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48054, 0, 3,
                                                                       47574, 21798, 47584, 5625,
                                                                       5634, 22086, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48084, 0, 3,
                                                                       47584, 21804, 47594, 5634,
                                                                       5643, 22104, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48114, 0, 3,
                                                                       47594, 21810, 47604, 5643,
                                                                       5652, 22122, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48144, 0, 3,
                                                                       47604, 21816, 47614, 5652,
                                                                       5661, 22140, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48174, 0, 3,
                                                                       47614, 21822, 47624, 5661,
                                                                       5670, 22158, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48204, 0, 3,
                                                                       47634, 21834, 47644, 5688,
                                                                       5697, 22176, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48234, 0, 3,
                                                                       47644, 21840, 47654, 5697,
                                                                       5706, 22194, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48264, 0, 3,
                                                                       47654, 21846, 47664, 5706,
                                                                       5715, 22212, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48294, 0, 3,
                                                                       47664, 21852, 47674, 5715,
                                                                       5724, 22230, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48324, 0, 3,
                                                                       47674, 21858, 47684, 5724,
                                                                       5733, 22248, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48354, 0, 3,
                                                                       47684, 21864, 47694, 5733,
                                                                       5742, 22266, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48384, 0, 3,
                                                                       47694, 21870, 47704, 5742,
                                                                       5751, 22284, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48414, 0, 3,
                                                                       47704, 21876, 47714, 5751,
                                                                       5760, 22302, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48444, 0, 3,
                                                                       47714, 21882, 47724, 5760,
                                                                       5769, 22320, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48474, 0, 3,
                                                                       47724, 21888, 47734, 5769,
                                                                       5778, 22338, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48504, 0, 3,
                                                                       47734, 21894, 47744, 5778,
                                                                       5787, 22356, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48534, 0, 3,
                                                                       47744, 21900, 47754, 5787,
                                                                       5796, 22374, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48564, 0, 3,
                                                                       47754, 21906, 47764, 5796,
                                                                       5805, 22392, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48594, 0, 3,
                                                                       47764, 21912, 47774, 5805,
                                                                       5814, 22410, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48624, 0, 3,
                                                                       47784, 21924, 47814, 5832,
                                                                       5850, 22428, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48684, 0, 3,
                                                                       47814, 21942, 47844, 5850,
                                                                       5868, 22464, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48744, 0, 3,
                                                                       47844, 21960, 47874, 5868,
                                                                       5886, 22500, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48804, 0, 3,
                                                                       47874, 21978, 47904, 5886,
                                                                       5904, 22536, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48864, 0, 3,
                                                                       47904, 21996, 47934, 5904,
                                                                       5922, 22572, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48924, 0, 3,
                                                                       47934, 22014, 47964, 5922,
                                                                       5940, 22608, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48984, 0, 3,
                                                                       47964, 22032, 47994, 5940,
                                                                       5958, 22644, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49044, 0, 3,
                                                                       47994, 22050, 48024, 5958,
                                                                       5976, 22680, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49104, 0, 3,
                                                                       48024, 22068, 48054, 5976,
                                                                       5994, 22716, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49164, 0, 3,
                                                                       48054, 22086, 48084, 5994,
                                                                       6012, 22752, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49224, 0, 3,
                                                                       48084, 22104, 48114, 6012,
                                                                       6030, 22788, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49284, 0, 3,
                                                                       48114, 22122, 48144, 6030,
                                                                       6048, 22824, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49344, 0, 3,
                                                                       48144, 22140, 48174, 6048,
                                                                       6066, 22860, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49404, 0, 3,
                                                                       48204, 22176, 48234, 6102,
                                                                       6120, 22896, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49464, 0, 3,
                                                                       48234, 22194, 48264, 6120,
                                                                       6138, 22932, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49524, 0, 3,
                                                                       48264, 22212, 48294, 6138,
                                                                       6156, 22968, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49584, 0, 3,
                                                                       48294, 22230, 48324, 6156,
                                                                       6174, 23004, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49644, 0, 3,
                                                                       48324, 22248, 48354, 6174,
                                                                       6192, 23040, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49704, 0, 3,
                                                                       48354, 22266, 48384, 6192,
                                                                       6210, 23076, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49764, 0, 3,
                                                                       48384, 22284, 48414, 6210,
                                                                       6228, 23112, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49824, 0, 3,
                                                                       48414, 22302, 48444, 6228,
                                                                       6246, 23148, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49884, 0, 3,
                                                                       48444, 22320, 48474, 6246,
                                                                       6264, 23184, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49944, 0, 3,
                                                                       48474, 22338, 48504, 6264,
                                                                       6282, 23220, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 50004, 0, 3,
                                                                       48504, 22356, 48534, 6282,
                                                                       6300, 23256, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 50064, 0, 3,
                                                                       48534, 22374, 48564, 6300,
                                                                       6318, 23292, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 50124, 0, 3,
                                                                       48564, 22392, 48594, 6318,
                                                                       6336, 23328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50184, 0, 3,
                                                                       48624, 22428, 48684, 6372,
                                                                       6402, 23364, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50284, 0, 3,
                                                                       48684, 22464, 48744, 6402,
                                                                       6432, 23424, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50384, 0, 3,
                                                                       48744, 22500, 48804, 6432,
                                                                       6462, 23484, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50484, 0, 3,
                                                                       48804, 22536, 48864, 6462,
                                                                       6492, 23544, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50584, 0, 3,
                                                                       48864, 22572, 48924, 6492,
                                                                       6522, 23604, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50684, 0, 3,
                                                                       48924, 22608, 48984, 6522,
                                                                       6552, 23664, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50784, 0, 3,
                                                                       48984, 22644, 49044, 6552,
                                                                       6582, 23724, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50884, 0, 3,
                                                                       49044, 22680, 49104, 6582,
                                                                       6612, 23784, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50984, 0, 3,
                                                                       49104, 22716, 49164, 6612,
                                                                       6642, 23844, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51084, 0, 3,
                                                                       49164, 22752, 49224, 6642,
                                                                       6672, 23904, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51184, 0, 3,
                                                                       49224, 22788, 49284, 6672,
                                                                       6702, 23964, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51284, 0, 3,
                                                                       49284, 22824, 49344, 6702,
                                                                       6732, 24024, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51384, 0, 3,
                                                                       49404, 22896, 49464, 6792,
                                                                       6822, 24084, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51484, 0, 3,
                                                                       49464, 22932, 49524, 6822,
                                                                       6852, 24144, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51584, 0, 3,
                                                                       49524, 22968, 49584, 6852,
                                                                       6882, 24204, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51684, 0, 3,
                                                                       49584, 23004, 49644, 6882,
                                                                       6912, 24264, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51784, 0, 3,
                                                                       49644, 23040, 49704, 6912,
                                                                       6942, 24324, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51884, 0, 3,
                                                                       49704, 23076, 49764, 6942,
                                                                       6972, 24384, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51984, 0, 3,
                                                                       49764, 23112, 49824, 6972,
                                                                       7002, 24444, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 52084, 0, 3,
                                                                       49824, 23148, 49884, 7002,
                                                                       7032, 24504, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 52184, 0, 3,
                                                                       49884, 23184, 49944, 7032,
                                                                       7062, 24564, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 52284, 0, 3,
                                                                       49944, 23220, 50004, 7062,
                                                                       7092, 24624, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 52384, 0, 3,
                                                                       50004, 23256, 50064, 7092,
                                                                       7122, 24684, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 52484, 0, 3,
                                                                       50064, 23292, 50124, 7122,
                                                                       7152, 24744, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 52584, 0, 3,
                                                                       50184, 23364, 50284, 7212,
                                                                       7257, 24804, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 52734, 0, 3,
                                                                       50284, 23424, 50384, 7257,
                                                                       7302, 24894, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 52884, 0, 3,
                                                                       50384, 23484, 50484, 7302,
                                                                       7347, 24984, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53034, 0, 3,
                                                                       50484, 23544, 50584, 7347,
                                                                       7392, 25074, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53184, 0, 3,
                                                                       50584, 23604, 50684, 7392,
                                                                       7437, 25164, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53334, 0, 3,
                                                                       50684, 23664, 50784, 7437,
                                                                       7482, 25254, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53484, 0, 3,
                                                                       50784, 23724, 50884, 7482,
                                                                       7527, 25344, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53634, 0, 3,
                                                                       50884, 23784, 50984, 7527,
                                                                       7572, 25434, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53784, 0, 3,
                                                                       50984, 23844, 51084, 7572,
                                                                       7617, 25524, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53934, 0, 3,
                                                                       51084, 23904, 51184, 7617,
                                                                       7662, 25614, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54084, 0, 3,
                                                                       51184, 23964, 51284, 7662,
                                                                       7707, 25704, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54234, 0, 3,
                                                                       51384, 24084, 51484, 7797,
                                                                       7842, 25794, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54384, 0, 3,
                                                                       51484, 24144, 51584, 7842,
                                                                       7887, 25884, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54534, 0, 3,
                                                                       51584, 24204, 51684, 7887,
                                                                       7932, 25974, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54684, 0, 3,
                                                                       51684, 24264, 51784, 7932,
                                                                       7977, 26064, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54834, 0, 3,
                                                                       51784, 24324, 51884, 7977,
                                                                       8022, 26154, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54984, 0, 3,
                                                                       51884, 24384, 51984, 8022,
                                                                       8067, 26244, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 55134, 0, 3,
                                                                       51984, 24444, 52084, 8067,
                                                                       8112, 26334, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 55284, 0, 3,
                                                                       52084, 24504, 52184, 8112,
                                                                       8157, 26424, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 55434, 0, 3,
                                                                       52184, 24564, 52284, 8157,
                                                                       8202, 26514, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 55584, 0, 3,
                                                                       52284, 24624, 52384, 8202,
                                                                       8247, 26604, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 55734, 0, 3,
                                                                       52384, 24684, 52484, 8247,
                                                                       8292, 26694, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 55884, 0, 3,
                                                                       52584, 24804, 52734, 8382,
                                                                       8445, 26784, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 56094, 0, 3,
                                                                       52734, 24894, 52884, 8445,
                                                                       8508, 26910, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 56304, 0, 3,
                                                                       52884, 24984, 53034, 8508,
                                                                       8571, 27036, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 56514, 0, 3,
                                                                       53034, 25074, 53184, 8571,
                                                                       8634, 27162, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 56724, 0, 3,
                                                                       53184, 25164, 53334, 8634,
                                                                       8697, 27288, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 56934, 0, 3,
                                                                       53334, 25254, 53484, 8697,
                                                                       8760, 27414, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 57144, 0, 3,
                                                                       53484, 25344, 53634, 8760,
                                                                       8823, 27540, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 57354, 0, 3,
                                                                       53634, 25434, 53784, 8823,
                                                                       8886, 27666, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 57564, 0, 3,
                                                                       53784, 25524, 53934, 8886,
                                                                       8949, 27792, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 57774, 0, 3,
                                                                       53934, 25614, 54084, 8949,
                                                                       9012, 27918, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 57984, 0, 3,
                                                                       54234, 25794, 54384, 9138,
                                                                       9201, 28044, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 58194, 0, 3,
                                                                       54384, 25884, 54534, 9201,
                                                                       9264, 28170, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 58404, 0, 3,
                                                                       54534, 25974, 54684, 9264,
                                                                       9327, 28296, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 58614, 0, 3,
                                                                       54684, 26064, 54834, 9327,
                                                                       9390, 28422, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 58824, 0, 3,
                                                                       54834, 26154, 54984, 9390,
                                                                       9453, 28548, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 59034, 0, 3,
                                                                       54984, 26244, 55134, 9453,
                                                                       9516, 28674, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 59244, 0, 3,
                                                                       55134, 26334, 55284, 9516,
                                                                       9579, 28800, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 59454, 0, 3,
                                                                       55284, 26424, 55434, 9579,
                                                                       9642, 28926, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 59664, 0, 3,
                                                                       55434, 26514, 55584, 9642,
                                                                       9705, 29052, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 59874, 0, 3,
                                                                       55584, 26604, 55734, 9705,
                                                                       9768, 29178, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 60084, 0, 3,
                                                                       55884, 26784, 56094, 9894,
                                                                       9978, 29304, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 60364, 0, 3,
                                                                       56094, 26910, 56304, 9978,
                                                                       10062, 29472, ncols,
                                                                       gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 60644, 0, 3,
                                                                       56304, 27036, 56514,
                                                                       10062, 10146, 29640,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 60924, 0, 3,
                                                                       56514, 27162, 56724,
                                                                       10146, 10230, 29808,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 61204, 0, 3,
                                                                       56724, 27288, 56934,
                                                                       10230, 10314, 29976,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 61484, 0, 3,
                                                                       56934, 27414, 57144,
                                                                       10314, 10398, 30144,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 61764, 0, 3,
                                                                       57144, 27540, 57354,
                                                                       10398, 10482, 30312,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 62044, 0, 3,
                                                                       57354, 27666, 57564,
                                                                       10482, 10566, 30480,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 62324, 0, 3,
                                                                       57564, 27792, 57774,
                                                                       10566, 10650, 30648,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 62604, 0, 3,
                                                                       57984, 28044, 58194,
                                                                       10818, 10902, 30816,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 62884, 0, 3,
                                                                       58194, 28170, 58404,
                                                                       10902, 10986, 30984,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 63164, 0, 3,
                                                                       58404, 28296, 58614,
                                                                       10986, 11070, 31152,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 63444, 0, 3,
                                                                       58614, 28422, 58824,
                                                                       11070, 11154, 31320,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 63724, 0, 3,
                                                                       58824, 28548, 59034,
                                                                       11154, 11238, 31488,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 64004, 0, 3,
                                                                       59034, 28674, 59244,
                                                                       11238, 11322, 31656,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 64284, 0, 3,
                                                                       59244, 28800, 59454,
                                                                       11322, 11406, 31824,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 64564, 0, 3,
                                                                       59454, 28926, 59664,
                                                                       11406, 11490, 31992,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 64844, 0, 3,
                                                                       59664, 29052, 59874,
                                                                       11490, 11574, 32160,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 65124, 0, 3,
                                                                       60084, 29304, 60364,
                                                                       11742, 11850, 32328,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 65484, 0, 3,
                                                                       60364, 29472, 60644,
                                                                       11850, 11958, 32544,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 65844, 0, 3,
                                                                       60644, 29640, 60924,
                                                                       11958, 12066, 32760,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 66204, 0, 3,
                                                                       60924, 29808, 61204,
                                                                       12066, 12174, 32976,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 66564, 0, 3,
                                                                       61204, 29976, 61484,
                                                                       12174, 12282, 33192,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 66924, 0, 3,
                                                                       61484, 30144, 61764,
                                                                       12282, 12390, 33408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 67284, 0, 3,
                                                                       61764, 30312, 62044,
                                                                       12390, 12498, 33624,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 67644, 0, 3,
                                                                       62044, 30480, 62324,
                                                                       12498, 12606, 33840,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 68004, 0, 3,
                                                                       62604, 30816, 62884,
                                                                       12822, 12930, 34056,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 68364, 0, 3,
                                                                       62884, 30984, 63164,
                                                                       12930, 13038, 34272,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 68724, 0, 3,
                                                                       63164, 31152, 63444,
                                                                       13038, 13146, 34488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 69084, 0, 3,
                                                                       63444, 31320, 63724,
                                                                       13146, 13254, 34704,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 69444, 0, 3,
                                                                       63724, 31488, 64004,
                                                                       13254, 13362, 34920,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 69804, 0, 3,
                                                                       64004, 31656, 64284,
                                                                       13362, 13470, 35136,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 70164, 0, 3,
                                                                       64284, 31824, 64564,
                                                                       13470, 13578, 35352,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 70524, 0, 3,
                                                                       64564, 31992, 64844,
                                                                       13578, 13686, 35568,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 70884, 0, 3,
                                                                       65124, 32328, 65484,
                                                                       13902, 14037, 35784,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 71334, 0, 3,
                                                                       65484, 32544, 65844,
                                                                       14037, 14172, 36054,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 71784, 0, 3,
                                                                       65844, 32760, 66204,
                                                                       14172, 14307, 36324,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 72234, 0, 3,
                                                                       66204, 32976, 66564,
                                                                       14307, 14442, 36594,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 72684, 0, 3,
                                                                       66564, 33192, 66924,
                                                                       14442, 14577, 36864,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 73134, 0, 3,
                                                                       66924, 33408, 67284,
                                                                       14577, 14712, 37134,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 73584, 0, 3,
                                                                       67284, 33624, 67644,
                                                                       14712, 14847, 37404,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 74034, 0, 3,
                                                                       68004, 34056, 68364,
                                                                       15117, 15252, 37674,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 74484, 0, 3,
                                                                       68364, 34272, 68724,
                                                                       15252, 15387, 37944,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 74934, 0, 3,
                                                                       68724, 34488, 69084,
                                                                       15387, 15522, 38214,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 75384, 0, 3,
                                                                       69084, 34704, 69444,
                                                                       15522, 15657, 38484,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 75834, 0, 3,
                                                                       69444, 34920, 69804,
                                                                       15657, 15792, 38754,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 76284, 0, 3,
                                                                       69804, 35136, 70164,
                                                                       15792, 15927, 39024,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 76734, 0, 3,
                                                                       70164, 35352, 70524,
                                                                       15927, 16062, 39294,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 77184, 0, 3,
                                                                       70884, 35784, 71334,
                                                                       16332, 16497, 39564,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 77734, 0, 3,
                                                                       71334, 36054, 71784,
                                                                       16497, 16662, 39894,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 78284, 0, 3,
                                                                       71784, 36324, 72234,
                                                                       16662, 16827, 40224,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 78834, 0, 3,
                                                                       72234, 36594, 72684,
                                                                       16827, 16992, 40554,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 79384, 0, 3,
                                                                       72684, 36864, 73134,
                                                                       16992, 17157, 40884,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 79934, 0, 3,
                                                                       73134, 37134, 73584,
                                                                       17157, 17322, 41214,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 80484, 0, 3,
                                                                       74034, 37674, 74484,
                                                                       17652, 17817, 41544,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 81034, 0, 3,
                                                                       74484, 37944, 74934,
                                                                       17817, 17982, 41874,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 81584, 0, 3,
                                                                       74934, 38214, 75384,
                                                                       17982, 18147, 42204,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 82134, 0, 3,
                                                                       75384, 38484, 75834,
                                                                       18147, 18312, 42534,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 82684, 0, 3,
                                                                       75834, 38754, 76284,
                                                                       18312, 18477, 42864,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 83234, 0, 3,
                                                                       76284, 39024, 76734,
                                                                       18477, 18642, 43194,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 83784, 0, 3,
                                                                       77184, 39564, 77734,
                                                                       18972, 19170, 43524,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 84444, 0, 3,
                                                                       77734, 39894, 78284,
                                                                       19170, 19368, 43920,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 85104, 0, 3,
                                                                       78284, 40224, 78834,
                                                                       19368, 19566, 44316,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 85764, 0, 3,
                                                                       78834, 40554, 79384,
                                                                       19566, 19764, 44712,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 86424, 0, 3,
                                                                       79384, 40884, 79934,
                                                                       19764, 19962, 45108,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 87084, 0, 3,
                                                                       80484, 41544, 81034,
                                                                       20358, 20556, 45504,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 87744, 0, 3,
                                                                       81034, 41874, 81584,
                                                                       20556, 20754, 45900,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 88404, 0, 3,
                                                                       81584, 42204, 82134,
                                                                       20754, 20952, 46296,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 89064, 0, 3,
                                                                       82134, 42534, 82684,
                                                                       20952, 21150, 46692,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 89724, 0, 3,
                                                                       82684, 42864, 83234,
                                                                       21150, 21348, 47088,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90384, 3, 21744,
                                                                       21750, 47504, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90399, 3, 21750,
                                                                       21756, 47514, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90414, 3, 21756,
                                                                       21762, 47524, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90429, 3, 21762,
                                                                       21768, 47534, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90444, 3, 21768,
                                                                       21774, 47544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90459, 3, 21774,
                                                                       21780, 47554, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90474, 3, 21780,
                                                                       21786, 47564, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90489, 3, 21786,
                                                                       21792, 47574, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90504, 3, 21792,
                                                                       21798, 47584, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90519, 3, 21798,
                                                                       21804, 47594, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90534, 3, 21804,
                                                                       21810, 47604, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90549, 3, 21810,
                                                                       21816, 47614, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90564, 3, 21816,
                                                                       21822, 47624, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90579, 3, 21834,
                                                                       21840, 47654, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90594, 3, 21840,
                                                                       21846, 47664, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90609, 3, 21846,
                                                                       21852, 47674, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90624, 3, 21852,
                                                                       21858, 47684, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90639, 3, 21858,
                                                                       21864, 47694, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90654, 3, 21864,
                                                                       21870, 47704, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90669, 3, 21870,
                                                                       21876, 47714, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90684, 3, 21876,
                                                                       21882, 47724, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90699, 3, 21882,
                                                                       21888, 47734, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90714, 3, 21888,
                                                                       21894, 47744, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90729, 3, 21894,
                                                                       21900, 47754, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90744, 3, 21900,
                                                                       21906, 47764, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90759, 3, 21906,
                                                                       21912, 47774, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90774, 0, 3,
                                                                       90384, 47504, 90399,
                                                                       21924, 21942, 47844,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90819, 0, 3,
                                                                       90399, 47514, 90414,
                                                                       21942, 21960, 47874,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90864, 0, 3,
                                                                       90414, 47524, 90429,
                                                                       21960, 21978, 47904,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90909, 0, 3,
                                                                       90429, 47534, 90444,
                                                                       21978, 21996, 47934,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90954, 0, 3,
                                                                       90444, 47544, 90459,
                                                                       21996, 22014, 47964,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90999, 0, 3,
                                                                       90459, 47554, 90474,
                                                                       22014, 22032, 47994,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91044, 0, 3,
                                                                       90474, 47564, 90489,
                                                                       22032, 22050, 48024,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91089, 0, 3,
                                                                       90489, 47574, 90504,
                                                                       22050, 22068, 48054,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91134, 0, 3,
                                                                       90504, 47584, 90519,
                                                                       22068, 22086, 48084,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91179, 0, 3,
                                                                       90519, 47594, 90534,
                                                                       22086, 22104, 48114,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91224, 0, 3,
                                                                       90534, 47604, 90549,
                                                                       22104, 22122, 48144,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91269, 0, 3,
                                                                       90549, 47614, 90564,
                                                                       22122, 22140, 48174,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91314, 0, 3,
                                                                       90579, 47654, 90594,
                                                                       22176, 22194, 48264,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91359, 0, 3,
                                                                       90594, 47664, 90609,
                                                                       22194, 22212, 48294,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91404, 0, 3,
                                                                       90609, 47674, 90624,
                                                                       22212, 22230, 48324,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91449, 0, 3,
                                                                       90624, 47684, 90639,
                                                                       22230, 22248, 48354,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91494, 0, 3,
                                                                       90639, 47694, 90654,
                                                                       22248, 22266, 48384,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91539, 0, 3,
                                                                       90654, 47704, 90669,
                                                                       22266, 22284, 48414,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91584, 0, 3,
                                                                       90669, 47714, 90684,
                                                                       22284, 22302, 48444,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91629, 0, 3,
                                                                       90684, 47724, 90699,
                                                                       22302, 22320, 48474,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91674, 0, 3,
                                                                       90699, 47734, 90714,
                                                                       22320, 22338, 48504,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91719, 0, 3,
                                                                       90714, 47744, 90729,
                                                                       22338, 22356, 48534,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91764, 0, 3,
                                                                       90729, 47754, 90744,
                                                                       22356, 22374, 48564,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 91809, 0, 3,
                                                                       90744, 47764, 90759,
                                                                       22374, 22392, 48594,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91854, 0, 3,
                                                                       90774, 47844, 90819,
                                                                       22428, 22464, 48744,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91944, 0, 3,
                                                                       90819, 47874, 90864,
                                                                       22464, 22500, 48804,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92034, 0, 3,
                                                                       90864, 47904, 90909,
                                                                       22500, 22536, 48864,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92124, 0, 3,
                                                                       90909, 47934, 90954,
                                                                       22536, 22572, 48924,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92214, 0, 3,
                                                                       90954, 47964, 90999,
                                                                       22572, 22608, 48984,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92304, 0, 3,
                                                                       90999, 47994, 91044,
                                                                       22608, 22644, 49044,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92394, 0, 3,
                                                                       91044, 48024, 91089,
                                                                       22644, 22680, 49104,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92484, 0, 3,
                                                                       91089, 48054, 91134,
                                                                       22680, 22716, 49164,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92574, 0, 3,
                                                                       91134, 48084, 91179,
                                                                       22716, 22752, 49224,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92664, 0, 3,
                                                                       91179, 48114, 91224,
                                                                       22752, 22788, 49284,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92754, 0, 3,
                                                                       91224, 48144, 91269,
                                                                       22788, 22824, 49344,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92844, 0, 3,
                                                                       91314, 48264, 91359,
                                                                       22896, 22932, 49524,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92934, 0, 3,
                                                                       91359, 48294, 91404,
                                                                       22932, 22968, 49584,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 93024, 0, 3,
                                                                       91404, 48324, 91449,
                                                                       22968, 23004, 49644,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 93114, 0, 3,
                                                                       91449, 48354, 91494,
                                                                       23004, 23040, 49704,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 93204, 0, 3,
                                                                       91494, 48384, 91539,
                                                                       23040, 23076, 49764,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 93294, 0, 3,
                                                                       91539, 48414, 91584,
                                                                       23076, 23112, 49824,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 93384, 0, 3,
                                                                       91584, 48444, 91629,
                                                                       23112, 23148, 49884,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 93474, 0, 3,
                                                                       91629, 48474, 91674,
                                                                       23148, 23184, 49944,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 93564, 0, 3,
                                                                       91674, 48504, 91719,
                                                                       23184, 23220, 50004,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 93654, 0, 3,
                                                                       91719, 48534, 91764,
                                                                       23220, 23256, 50064,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 93744, 0, 3,
                                                                       91764, 48564, 91809,
                                                                       23256, 23292, 50124,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 93834, 0, 3,
                                                                       91854, 48744, 91944,
                                                                       23364, 23424, 50384,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 93984, 0, 3,
                                                                       91944, 48804, 92034,
                                                                       23424, 23484, 50484,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94134, 0, 3,
                                                                       92034, 48864, 92124,
                                                                       23484, 23544, 50584,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94284, 0, 3,
                                                                       92124, 48924, 92214,
                                                                       23544, 23604, 50684,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94434, 0, 3,
                                                                       92214, 48984, 92304,
                                                                       23604, 23664, 50784,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94584, 0, 3,
                                                                       92304, 49044, 92394,
                                                                       23664, 23724, 50884,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94734, 0, 3,
                                                                       92394, 49104, 92484,
                                                                       23724, 23784, 50984,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94884, 0, 3,
                                                                       92484, 49164, 92574,
                                                                       23784, 23844, 51084,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 95034, 0, 3,
                                                                       92574, 49224, 92664,
                                                                       23844, 23904, 51184,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 95184, 0, 3,
                                                                       92664, 49284, 92754,
                                                                       23904, 23964, 51284,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 95334, 0, 3,
                                                                       92844, 49524, 92934,
                                                                       24084, 24144, 51584,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 95484, 0, 3,
                                                                       92934, 49584, 93024,
                                                                       24144, 24204, 51684,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 95634, 0, 3,
                                                                       93024, 49644, 93114,
                                                                       24204, 24264, 51784,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 95784, 0, 3,
                                                                       93114, 49704, 93204,
                                                                       24264, 24324, 51884,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 95934, 0, 3,
                                                                       93204, 49764, 93294,
                                                                       24324, 24384, 51984,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 96084, 0, 3,
                                                                       93294, 49824, 93384,
                                                                       24384, 24444, 52084,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 96234, 0, 3,
                                                                       93384, 49884, 93474,
                                                                       24444, 24504, 52184,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 96384, 0, 3,
                                                                       93474, 49944, 93564,
                                                                       24504, 24564, 52284,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 96534, 0, 3,
                                                                       93564, 50004, 93654,
                                                                       24564, 24624, 52384,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 96684, 0, 3,
                                                                       93654, 50064, 93744,
                                                                       24624, 24684, 52484,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 96834, 0, 3,
                                                                       93834, 50384, 93984,
                                                                       24804, 24894, 52884,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 97059, 0, 3,
                                                                       93984, 50484, 94134,
                                                                       24894, 24984, 53034,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 97284, 0, 3,
                                                                       94134, 50584, 94284,
                                                                       24984, 25074, 53184,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 97509, 0, 3,
                                                                       94284, 50684, 94434,
                                                                       25074, 25164, 53334,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 97734, 0, 3,
                                                                       94434, 50784, 94584,
                                                                       25164, 25254, 53484,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 97959, 0, 3,
                                                                       94584, 50884, 94734,
                                                                       25254, 25344, 53634,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 98184, 0, 3,
                                                                       94734, 50984, 94884,
                                                                       25344, 25434, 53784,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 98409, 0, 3,
                                                                       94884, 51084, 95034,
                                                                       25434, 25524, 53934,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 98634, 0, 3,
                                                                       95034, 51184, 95184,
                                                                       25524, 25614, 54084,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 98859, 0, 3,
                                                                       95334, 51584, 95484,
                                                                       25794, 25884, 54534,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 99084, 0, 3,
                                                                       95484, 51684, 95634,
                                                                       25884, 25974, 54684,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 99309, 0, 3,
                                                                       95634, 51784, 95784,
                                                                       25974, 26064, 54834,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 99534, 0, 3,
                                                                       95784, 51884, 95934,
                                                                       26064, 26154, 54984,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 99759, 0, 3,
                                                                       95934, 51984, 96084,
                                                                       26154, 26244, 55134,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 99984, 0, 3,
                                                                       96084, 52084, 96234,
                                                                       26244, 26334, 55284,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 100209, 0, 3,
                                                                       96234, 52184, 96384,
                                                                       26334, 26424, 55434,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 100434, 0, 3,
                                                                       96384, 52284, 96534,
                                                                       26424, 26514, 55584,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 100659, 0, 3,
                                                                       96534, 52384, 96684,
                                                                       26514, 26604, 55734,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 100884, 0, 3,
                                                                       96834, 52884, 97059,
                                                                       26784, 26910, 56304,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 101199, 0, 3,
                                                                       97059, 53034, 97284,
                                                                       26910, 27036, 56514,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 101514, 0, 3,
                                                                       97284, 53184, 97509,
                                                                       27036, 27162, 56724,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 101829, 0, 3,
                                                                       97509, 53334, 97734,
                                                                       27162, 27288, 56934,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 102144, 0, 3,
                                                                       97734, 53484, 97959,
                                                                       27288, 27414, 57144,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 102459, 0, 3,
                                                                       97959, 53634, 98184,
                                                                       27414, 27540, 57354,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 102774, 0, 3,
                                                                       98184, 53784, 98409,
                                                                       27540, 27666, 57564,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 103089, 0, 3,
                                                                       98409, 53934, 98634,
                                                                       27666, 27792, 57774,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 103404, 0, 3,
                                                                       98859, 54534, 99084,
                                                                       28044, 28170, 58404,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 103719, 0, 3,
                                                                       99084, 54684, 99309,
                                                                       28170, 28296, 58614,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 104034, 0, 3,
                                                                       99309, 54834, 99534,
                                                                       28296, 28422, 58824,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 104349, 0, 3,
                                                                       99534, 54984, 99759,
                                                                       28422, 28548, 59034,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 104664, 0, 3,
                                                                       99759, 55134, 99984,
                                                                       28548, 28674, 59244,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 104979, 0, 3,
                                                                       99984, 55284, 100209,
                                                                       28674, 28800, 59454,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 105294, 0, 3,
                                                                       100209, 55434, 100434,
                                                                       28800, 28926, 59664,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 105609, 0, 3,
                                                                       100434, 55584, 100659,
                                                                       28926, 29052, 59874,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 105924, 0, 3,
                                                                       100884, 56304, 101199,
                                                                       29304, 29472, 60644,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 106344, 0, 3,
                                                                       101199, 56514, 101514,
                                                                       29472, 29640, 60924,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 106764, 0, 3,
                                                                       101514, 56724, 101829,
                                                                       29640, 29808, 61204,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 107184, 0, 3,
                                                                       101829, 56934, 102144,
                                                                       29808, 29976, 61484,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 107604, 0, 3,
                                                                       102144, 57144, 102459,
                                                                       29976, 30144, 61764,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 108024, 0, 3,
                                                                       102459, 57354, 102774,
                                                                       30144, 30312, 62044,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 108444, 0, 3,
                                                                       102774, 57564, 103089,
                                                                       30312, 30480, 62324,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 108864, 0, 3,
                                                                       103404, 58404, 103719,
                                                                       30816, 30984, 63164,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 109284, 0, 3,
                                                                       103719, 58614, 104034,
                                                                       30984, 31152, 63444,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 109704, 0, 3,
                                                                       104034, 58824, 104349,
                                                                       31152, 31320, 63724,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 110124, 0, 3,
                                                                       104349, 59034, 104664,
                                                                       31320, 31488, 64004,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 110544, 0, 3,
                                                                       104664, 59244, 104979,
                                                                       31488, 31656, 64284,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 110964, 0, 3,
                                                                       104979, 59454, 105294,
                                                                       31656, 31824, 64564,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 111384, 0, 3,
                                                                       105294, 59664, 105609,
                                                                       31824, 31992, 64844,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 111804, 0, 3,
                                                                       105924, 60644, 106344,
                                                                       32328, 32544, 65844,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 112344, 0, 3,
                                                                       106344, 60924, 106764,
                                                                       32544, 32760, 66204,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 112884, 0, 3,
                                                                       106764, 61204, 107184,
                                                                       32760, 32976, 66564,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 113424, 0, 3,
                                                                       107184, 61484, 107604,
                                                                       32976, 33192, 66924,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 113964, 0, 3,
                                                                       107604, 61764, 108024,
                                                                       33192, 33408, 67284,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 114504, 0, 3,
                                                                       108024, 62044, 108444,
                                                                       33408, 33624, 67644,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 115044, 0, 3,
                                                                       108864, 63164, 109284,
                                                                       34056, 34272, 68724,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 115584, 0, 3,
                                                                       109284, 63444, 109704,
                                                                       34272, 34488, 69084,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 116124, 0, 3,
                                                                       109704, 63724, 110124,
                                                                       34488, 34704, 69444,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 116664, 0, 3,
                                                                       110124, 64004, 110544,
                                                                       34704, 34920, 69804,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 117204, 0, 3,
                                                                       110544, 64284, 110964,
                                                                       34920, 35136, 70164,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 117744, 0, 3,
                                                                       110964, 64564, 111384,
                                                                       35136, 35352, 70524,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 118284, 0, 3,
                                                                       111804, 65844, 112344,
                                                                       35784, 36054, 71784,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 118959, 0, 3,
                                                                       112344, 66204, 112884,
                                                                       36054, 36324, 72234,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 119634, 0, 3,
                                                                       112884, 66564, 113424,
                                                                       36324, 36594, 72684,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 120309, 0, 3,
                                                                       113424, 66924, 113964,
                                                                       36594, 36864, 73134,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 120984, 0, 3,
                                                                       113964, 67284, 114504,
                                                                       36864, 37134, 73584,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 121659, 0, 3,
                                                                       115044, 68724, 115584,
                                                                       37674, 37944, 74934,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 122334, 0, 3,
                                                                       115584, 69084, 116124,
                                                                       37944, 38214, 75384,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 123009, 0, 3,
                                                                       116124, 69444, 116664,
                                                                       38214, 38484, 75834,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 123684, 0, 3,
                                                                       116664, 69804, 117204,
                                                                       38484, 38754, 76284,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 124359, 0, 3,
                                                                       117204, 70164, 117744,
                                                                       38754, 39024, 76734,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 125034, 0, 3,
                                                                       118284, 71784, 118959,
                                                                       39564, 39894, 78284,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 125859, 0, 3,
                                                                       118959, 72234, 119634,
                                                                       39894, 40224, 78834,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 126684, 0, 3,
                                                                       119634, 72684, 120309,
                                                                       40224, 40554, 79384,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 127509, 0, 3,
                                                                       120309, 73134, 120984,
                                                                       40554, 40884, 79934,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 128334, 0, 3,
                                                                       121659, 74934, 122334,
                                                                       41544, 41874, 81584,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 129159, 0, 3,
                                                                       122334, 75384, 123009,
                                                                       41874, 42204, 82134,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 129984, 0, 3,
                                                                       123009, 75834, 123684,
                                                                       42204, 42534, 82684,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 130809, 0, 3,
                                                                       123684, 76284, 124359,
                                                                       42534, 42864, 83234,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 131634, 0, 3,
                                                                       125034, 78284, 125859,
                                                                       43524, 43920, 85104,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 132624, 0, 3,
                                                                       125859, 78834, 126684,
                                                                       43920, 44316, 85764,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 133614, 0, 3,
                                                                       126684, 79384, 127509,
                                                                       44316, 44712, 86424,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 134604, 0, 3,
                                                                       128334, 81584, 129159,
                                                                       45504, 45900, 88404,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 135594, 0, 3,
                                                                       129159, 82134, 129984,
                                                                       45900, 46296, 89064,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 136584, 0, 3,
                                                                       129984, 82684, 130809,
                                                                       46296, 46692, 89724,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137574, 3, 47484,
                                                                       47494, 90384, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137595, 3, 47494,
                                                                       47504, 90399, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137616, 3, 47504,
                                                                       47514, 90414, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137637, 3, 47514,
                                                                       47524, 90429, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137658, 3, 47524,
                                                                       47534, 90444, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137679, 3, 47534,
                                                                       47544, 90459, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137700, 3, 47544,
                                                                       47554, 90474, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137721, 3, 47554,
                                                                       47564, 90489, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137742, 3, 47564,
                                                                       47574, 90504, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137763, 3, 47574,
                                                                       47584, 90519, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137784, 3, 47584,
                                                                       47594, 90534, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137805, 3, 47594,
                                                                       47604, 90549, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137826, 3, 47604,
                                                                       47614, 90564, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137847, 3, 47634,
                                                                       47644, 90579, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137868, 3, 47644,
                                                                       47654, 90594, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137889, 3, 47654,
                                                                       47664, 90609, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137910, 3, 47664,
                                                                       47674, 90624, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137931, 3, 47674,
                                                                       47684, 90639, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137952, 3, 47684,
                                                                       47694, 90654, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137973, 3, 47694,
                                                                       47704, 90669, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137994, 3, 47704,
                                                                       47714, 90684, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 138015, 3, 47714,
                                                                       47724, 90699, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 138036, 3, 47724,
                                                                       47734, 90714, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 138057, 3, 47734,
                                                                       47744, 90729, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 138078, 3, 47744,
                                                                       47754, 90744, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 138099, 3, 47754,
                                                                       47764, 90759, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138120, 0, 3,
                                                                       137574, 90384, 137595,
                                                                       47784, 47814, 90774,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138183, 0, 3,
                                                                       137595, 90399, 137616,
                                                                       47814, 47844, 90819,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138246, 0, 3,
                                                                       137616, 90414, 137637,
                                                                       47844, 47874, 90864,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138309, 0, 3,
                                                                       137637, 90429, 137658,
                                                                       47874, 47904, 90909,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138372, 0, 3,
                                                                       137658, 90444, 137679,
                                                                       47904, 47934, 90954,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138435, 0, 3,
                                                                       137679, 90459, 137700,
                                                                       47934, 47964, 90999,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138498, 0, 3,
                                                                       137700, 90474, 137721,
                                                                       47964, 47994, 91044,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138561, 0, 3,
                                                                       137721, 90489, 137742,
                                                                       47994, 48024, 91089,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138624, 0, 3,
                                                                       137742, 90504, 137763,
                                                                       48024, 48054, 91134,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138687, 0, 3,
                                                                       137763, 90519, 137784,
                                                                       48054, 48084, 91179,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138750, 0, 3,
                                                                       137784, 90534, 137805,
                                                                       48084, 48114, 91224,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138813, 0, 3,
                                                                       137805, 90549, 137826,
                                                                       48114, 48144, 91269,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138876, 0, 3,
                                                                       137847, 90579, 137868,
                                                                       48204, 48234, 91314,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138939, 0, 3,
                                                                       137868, 90594, 137889,
                                                                       48234, 48264, 91359,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 139002, 0, 3,
                                                                       137889, 90609, 137910,
                                                                       48264, 48294, 91404,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 139065, 0, 3,
                                                                       137910, 90624, 137931,
                                                                       48294, 48324, 91449,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 139128, 0, 3,
                                                                       137931, 90639, 137952,
                                                                       48324, 48354, 91494,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 139191, 0, 3,
                                                                       137952, 90654, 137973,
                                                                       48354, 48384, 91539,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 139254, 0, 3,
                                                                       137973, 90669, 137994,
                                                                       48384, 48414, 91584,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 139317, 0, 3,
                                                                       137994, 90684, 138015,
                                                                       48414, 48444, 91629,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 139380, 0, 3,
                                                                       138015, 90699, 138036,
                                                                       48444, 48474, 91674,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 139443, 0, 3,
                                                                       138036, 90714, 138057,
                                                                       48474, 48504, 91719,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 139506, 0, 3,
                                                                       138057, 90729, 138078,
                                                                       48504, 48534, 91764,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 139569, 0, 3,
                                                                       138078, 90744, 138099,
                                                                       48534, 48564, 91809,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139632, 0, 3,
                                                                       138120, 90774, 138183,
                                                                       48624, 48684, 91854,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139758, 0, 3,
                                                                       138183, 90819, 138246,
                                                                       48684, 48744, 91944,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139884, 0, 3,
                                                                       138246, 90864, 138309,
                                                                       48744, 48804, 92034,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140010, 0, 3,
                                                                       138309, 90909, 138372,
                                                                       48804, 48864, 92124,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140136, 0, 3,
                                                                       138372, 90954, 138435,
                                                                       48864, 48924, 92214,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140262, 0, 3,
                                                                       138435, 90999, 138498,
                                                                       48924, 48984, 92304,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140388, 0, 3,
                                                                       138498, 91044, 138561,
                                                                       48984, 49044, 92394,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140514, 0, 3,
                                                                       138561, 91089, 138624,
                                                                       49044, 49104, 92484,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140640, 0, 3,
                                                                       138624, 91134, 138687,
                                                                       49104, 49164, 92574,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140766, 0, 3,
                                                                       138687, 91179, 138750,
                                                                       49164, 49224, 92664,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140892, 0, 3,
                                                                       138750, 91224, 138813,
                                                                       49224, 49284, 92754,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 141018, 0, 3,
                                                                       138876, 91314, 138939,
                                                                       49404, 49464, 92844,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 141144, 0, 3,
                                                                       138939, 91359, 139002,
                                                                       49464, 49524, 92934,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 141270, 0, 3,
                                                                       139002, 91404, 139065,
                                                                       49524, 49584, 93024,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 141396, 0, 3,
                                                                       139065, 91449, 139128,
                                                                       49584, 49644, 93114,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 141522, 0, 3,
                                                                       139128, 91494, 139191,
                                                                       49644, 49704, 93204,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 141648, 0, 3,
                                                                       139191, 91539, 139254,
                                                                       49704, 49764, 93294,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 141774, 0, 3,
                                                                       139254, 91584, 139317,
                                                                       49764, 49824, 93384,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 141900, 0, 3,
                                                                       139317, 91629, 139380,
                                                                       49824, 49884, 93474,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 142026, 0, 3,
                                                                       139380, 91674, 139443,
                                                                       49884, 49944, 93564,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 142152, 0, 3,
                                                                       139443, 91719, 139506,
                                                                       49944, 50004, 93654,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 142278, 0, 3,
                                                                       139506, 91764, 139569,
                                                                       50004, 50064, 93744,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 142404, 0, 3,
                                                                       139632, 91854, 139758,
                                                                       50184, 50284, 93834,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 142614, 0, 3,
                                                                       139758, 91944, 139884,
                                                                       50284, 50384, 93984,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 142824, 0, 3,
                                                                       139884, 92034, 140010,
                                                                       50384, 50484, 94134,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 143034, 0, 3,
                                                                       140010, 92124, 140136,
                                                                       50484, 50584, 94284,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 143244, 0, 3,
                                                                       140136, 92214, 140262,
                                                                       50584, 50684, 94434,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 143454, 0, 3,
                                                                       140262, 92304, 140388,
                                                                       50684, 50784, 94584,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 143664, 0, 3,
                                                                       140388, 92394, 140514,
                                                                       50784, 50884, 94734,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 143874, 0, 3,
                                                                       140514, 92484, 140640,
                                                                       50884, 50984, 94884,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 144084, 0, 3,
                                                                       140640, 92574, 140766,
                                                                       50984, 51084, 95034,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 144294, 0, 3,
                                                                       140766, 92664, 140892,
                                                                       51084, 51184, 95184,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 144504, 0, 3,
                                                                       141018, 92844, 141144,
                                                                       51384, 51484, 95334,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 144714, 0, 3,
                                                                       141144, 92934, 141270,
                                                                       51484, 51584, 95484,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 144924, 0, 3,
                                                                       141270, 93024, 141396,
                                                                       51584, 51684, 95634,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 145134, 0, 3,
                                                                       141396, 93114, 141522,
                                                                       51684, 51784, 95784,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 145344, 0, 3,
                                                                       141522, 93204, 141648,
                                                                       51784, 51884, 95934,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 145554, 0, 3,
                                                                       141648, 93294, 141774,
                                                                       51884, 51984, 96084,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 145764, 0, 3,
                                                                       141774, 93384, 141900,
                                                                       51984, 52084, 96234,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 145974, 0, 3,
                                                                       141900, 93474, 142026,
                                                                       52084, 52184, 96384,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 146184, 0, 3,
                                                                       142026, 93564, 142152,
                                                                       52184, 52284, 96534,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 146394, 0, 3,
                                                                       142152, 93654, 142278,
                                                                       52284, 52384, 96684,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 146604, 0, 3,
                                                                       142404, 93834, 142614,
                                                                       52584, 52734, 96834,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 146919, 0, 3,
                                                                       142614, 93984, 142824,
                                                                       52734, 52884, 97059,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 147234, 0, 3,
                                                                       142824, 94134, 143034,
                                                                       52884, 53034, 97284,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 147549, 0, 3,
                                                                       143034, 94284, 143244,
                                                                       53034, 53184, 97509,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 147864, 0, 3,
                                                                       143244, 94434, 143454,
                                                                       53184, 53334, 97734,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 148179, 0, 3,
                                                                       143454, 94584, 143664,
                                                                       53334, 53484, 97959,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 148494, 0, 3,
                                                                       143664, 94734, 143874,
                                                                       53484, 53634, 98184,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 148809, 0, 3,
                                                                       143874, 94884, 144084,
                                                                       53634, 53784, 98409,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 149124, 0, 3,
                                                                       144084, 95034, 144294,
                                                                       53784, 53934, 98634,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 149439, 0, 3,
                                                                       144504, 95334, 144714,
                                                                       54234, 54384, 98859,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 149754, 0, 3,
                                                                       144714, 95484, 144924,
                                                                       54384, 54534, 99084,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 150069, 0, 3,
                                                                       144924, 95634, 145134,
                                                                       54534, 54684, 99309,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 150384, 0, 3,
                                                                       145134, 95784, 145344,
                                                                       54684, 54834, 99534,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 150699, 0, 3,
                                                                       145344, 95934, 145554,
                                                                       54834, 54984, 99759,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 151014, 0, 3,
                                                                       145554, 96084, 145764,
                                                                       54984, 55134, 99984,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 151329, 0, 3,
                                                                       145764, 96234, 145974,
                                                                       55134, 55284, 100209,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 151644, 0, 3,
                                                                       145974, 96384, 146184,
                                                                       55284, 55434, 100434,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 151959, 0, 3,
                                                                       146184, 96534, 146394,
                                                                       55434, 55584, 100659,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 152274, 0, 3,
                                                                       146604, 96834, 146919,
                                                                       55884, 56094, 100884,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 152715, 0, 3,
                                                                       146919, 97059, 147234,
                                                                       56094, 56304, 101199,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 153156, 0, 3,
                                                                       147234, 97284, 147549,
                                                                       56304, 56514, 101514,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 153597, 0, 3,
                                                                       147549, 97509, 147864,
                                                                       56514, 56724, 101829,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 154038, 0, 3,
                                                                       147864, 97734, 148179,
                                                                       56724, 56934, 102144,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 154479, 0, 3,
                                                                       148179, 97959, 148494,
                                                                       56934, 57144, 102459,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 154920, 0, 3,
                                                                       148494, 98184, 148809,
                                                                       57144, 57354, 102774,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 155361, 0, 3,
                                                                       148809, 98409, 149124,
                                                                       57354, 57564, 103089,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 155802, 0, 3,
                                                                       149439, 98859, 149754,
                                                                       57984, 58194, 103404,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 156243, 0, 3,
                                                                       149754, 99084, 150069,
                                                                       58194, 58404, 103719,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 156684, 0, 3,
                                                                       150069, 99309, 150384,
                                                                       58404, 58614, 104034,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 157125, 0, 3,
                                                                       150384, 99534, 150699,
                                                                       58614, 58824, 104349,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 157566, 0, 3,
                                                                       150699, 99759, 151014,
                                                                       58824, 59034, 104664,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 158007, 0, 3,
                                                                       151014, 99984, 151329,
                                                                       59034, 59244, 104979,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 158448, 0, 3,
                                                                       151329, 100209, 151644,
                                                                       59244, 59454, 105294,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 158889, 0, 3,
                                                                       151644, 100434, 151959,
                                                                       59454, 59664, 105609,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 159330, 0, 3,
                                                                       152274, 100884, 152715,
                                                                       60084, 60364, 105924,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 159918, 0, 3,
                                                                       152715, 101199, 153156,
                                                                       60364, 60644, 106344,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 160506, 0, 3,
                                                                       153156, 101514, 153597,
                                                                       60644, 60924, 106764,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 161094, 0, 3,
                                                                       153597, 101829, 154038,
                                                                       60924, 61204, 107184,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 161682, 0, 3,
                                                                       154038, 102144, 154479,
                                                                       61204, 61484, 107604,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 162270, 0, 3,
                                                                       154479, 102459, 154920,
                                                                       61484, 61764, 108024,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 162858, 0, 3,
                                                                       154920, 102774, 155361,
                                                                       61764, 62044, 108444,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 163446, 0, 3,
                                                                       155802, 103404, 156243,
                                                                       62604, 62884, 108864,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 164034, 0, 3,
                                                                       156243, 103719, 156684,
                                                                       62884, 63164, 109284,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 164622, 0, 3,
                                                                       156684, 104034, 157125,
                                                                       63164, 63444, 109704,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 165210, 0, 3,
                                                                       157125, 104349, 157566,
                                                                       63444, 63724, 110124,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 165798, 0, 3,
                                                                       157566, 104664, 158007,
                                                                       63724, 64004, 110544,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 166386, 0, 3,
                                                                       158007, 104979, 158448,
                                                                       64004, 64284, 110964,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 166974, 0, 3,
                                                                       158448, 105294, 158889,
                                                                       64284, 64564, 111384,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 167562, 0, 3,
                                                                       159330, 105924, 159918,
                                                                       65124, 65484, 111804,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 168318, 0, 3,
                                                                       159918, 106344, 160506,
                                                                       65484, 65844, 112344,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 169074, 0, 3,
                                                                       160506, 106764, 161094,
                                                                       65844, 66204, 112884,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 169830, 0, 3,
                                                                       161094, 107184, 161682,
                                                                       66204, 66564, 113424,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 170586, 0, 3,
                                                                       161682, 107604, 162270,
                                                                       66564, 66924, 113964,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 171342, 0, 3,
                                                                       162270, 108024, 162858,
                                                                       66924, 67284, 114504,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 172098, 0, 3,
                                                                       163446, 108864, 164034,
                                                                       68004, 68364, 115044,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 172854, 0, 3,
                                                                       164034, 109284, 164622,
                                                                       68364, 68724, 115584,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 173610, 0, 3,
                                                                       164622, 109704, 165210,
                                                                       68724, 69084, 116124,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 174366, 0, 3,
                                                                       165210, 110124, 165798,
                                                                       69084, 69444, 116664,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 175122, 0, 3,
                                                                       165798, 110544, 166386,
                                                                       69444, 69804, 117204,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 175878, 0, 3,
                                                                       166386, 110964, 166974,
                                                                       69804, 70164, 117744,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 176634, 0, 3,
                                                                       167562, 111804, 168318,
                                                                       70884, 71334, 118284,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 177579, 0, 3,
                                                                       168318, 112344, 169074,
                                                                       71334, 71784, 118959,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 178524, 0, 3,
                                                                       169074, 112884, 169830,
                                                                       71784, 72234, 119634,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 179469, 0, 3,
                                                                       169830, 113424, 170586,
                                                                       72234, 72684, 120309,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 180414, 0, 3,
                                                                       170586, 113964, 171342,
                                                                       72684, 73134, 120984,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 181359, 0, 3,
                                                                       172098, 115044, 172854,
                                                                       74034, 74484, 121659,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 182304, 0, 3,
                                                                       172854, 115584, 173610,
                                                                       74484, 74934, 122334,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 183249, 0, 3,
                                                                       173610, 116124, 174366,
                                                                       74934, 75384, 123009,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 184194, 0, 3,
                                                                       174366, 116664, 175122,
                                                                       75384, 75834, 123684,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 185139, 0, 3,
                                                                       175122, 117204, 175878,
                                                                       75834, 76284, 124359,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 186084, 0, 3,
                                                                       176634, 118284, 177579,
                                                                       77184, 77734, 125034,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 187239, 0, 3,
                                                                       177579, 118959, 178524,
                                                                       77734, 78284, 125859,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 188394, 0, 3,
                                                                       178524, 119634, 179469,
                                                                       78284, 78834, 126684,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 189549, 0, 3,
                                                                       179469, 120309, 180414,
                                                                       78834, 79384, 127509,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 190704, 0, 3,
                                                                       181359, 121659, 182304,
                                                                       80484, 81034, 128334,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 191859, 0, 3,
                                                                       182304, 122334, 183249,
                                                                       81034, 81584, 129159,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 193014, 0, 3,
                                                                       183249, 123009, 184194,
                                                                       81584, 82134, 129984,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 194169, 0, 3,
                                                                       184194, 123684, 185139,
                                                                       82134, 82684, 130809,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 195324, 0, 3,
                                                                       186084, 125034, 187239,
                                                                       83784, 84444, 131634,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 196710, 0, 3,
                                                                       187239, 125859, 188394,
                                                                       84444, 85104, 132624,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 198096, 0, 3,
                                                                       188394, 126684, 189549,
                                                                       85104, 85764, 133614,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 199482, 0, 3,
                                                                       190704, 128334, 191859,
                                                                       87084, 87744, 134604,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 200868, 0, 3,
                                                                       191859, 129159, 193014,
                                                                       87744, 88404, 135594,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 202254, 0, 3,
                                                                       193014, 129984, 194169,
                                                                       88404, 89064, 136584,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203640, 3, 90384,
                                                                       90399, 137616, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203668, 3, 90399,
                                                                       90414, 137637, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203696, 3, 90414,
                                                                       90429, 137658, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203724, 3, 90429,
                                                                       90444, 137679, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203752, 3, 90444,
                                                                       90459, 137700, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203780, 3, 90459,
                                                                       90474, 137721, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203808, 3, 90474,
                                                                       90489, 137742, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203836, 3, 90489,
                                                                       90504, 137763, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203864, 3, 90504,
                                                                       90519, 137784, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203892, 3, 90519,
                                                                       90534, 137805, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203920, 3, 90534,
                                                                       90549, 137826, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203948, 3, 90579,
                                                                       90594, 137889, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203976, 3, 90594,
                                                                       90609, 137910, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204004, 3, 90609,
                                                                       90624, 137931, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204032, 3, 90624,
                                                                       90639, 137952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204060, 3, 90639,
                                                                       90654, 137973, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204088, 3, 90654,
                                                                       90669, 137994, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204116, 3, 90669,
                                                                       90684, 138015, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204144, 3, 90684,
                                                                       90699, 138036, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204172, 3, 90699,
                                                                       90714, 138057, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204200, 3, 90714,
                                                                       90729, 138078, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204228, 3, 90729,
                                                                       90744, 138099, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 204256, 0, 3,
                                                                       203640, 137616, 203668,
                                                                       90774, 90819, 138246,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 204340, 0, 3,
                                                                       203668, 137637, 203696,
                                                                       90819, 90864, 138309,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 204424, 0, 3,
                                                                       203696, 137658, 203724,
                                                                       90864, 90909, 138372,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 204508, 0, 3,
                                                                       203724, 137679, 203752,
                                                                       90909, 90954, 138435,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 204592, 0, 3,
                                                                       203752, 137700, 203780,
                                                                       90954, 90999, 138498,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 204676, 0, 3,
                                                                       203780, 137721, 203808,
                                                                       90999, 91044, 138561,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 204760, 0, 3,
                                                                       203808, 137742, 203836,
                                                                       91044, 91089, 138624,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 204844, 0, 3,
                                                                       203836, 137763, 203864,
                                                                       91089, 91134, 138687,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 204928, 0, 3,
                                                                       203864, 137784, 203892,
                                                                       91134, 91179, 138750,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 205012, 0, 3,
                                                                       203892, 137805, 203920,
                                                                       91179, 91224, 138813,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 205096, 0, 3,
                                                                       203948, 137889, 203976,
                                                                       91314, 91359, 139002,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 205180, 0, 3,
                                                                       203976, 137910, 204004,
                                                                       91359, 91404, 139065,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 205264, 0, 3,
                                                                       204004, 137931, 204032,
                                                                       91404, 91449, 139128,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 205348, 0, 3,
                                                                       204032, 137952, 204060,
                                                                       91449, 91494, 139191,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 205432, 0, 3,
                                                                       204060, 137973, 204088,
                                                                       91494, 91539, 139254,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 205516, 0, 3,
                                                                       204088, 137994, 204116,
                                                                       91539, 91584, 139317,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 205600, 0, 3,
                                                                       204116, 138015, 204144,
                                                                       91584, 91629, 139380,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 205684, 0, 3,
                                                                       204144, 138036, 204172,
                                                                       91629, 91674, 139443,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 205768, 0, 3,
                                                                       204172, 138057, 204200,
                                                                       91674, 91719, 139506,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 205852, 0, 3,
                                                                       204200, 138078, 204228,
                                                                       91719, 91764, 139569,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 205936, 0, 3,
                                                                       204256, 138246, 204340,
                                                                       91854, 91944, 139884,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 206104, 0, 3,
                                                                       204340, 138309, 204424,
                                                                       91944, 92034, 140010,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 206272, 0, 3,
                                                                       204424, 138372, 204508,
                                                                       92034, 92124, 140136,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 206440, 0, 3,
                                                                       204508, 138435, 204592,
                                                                       92124, 92214, 140262,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 206608, 0, 3,
                                                                       204592, 138498, 204676,
                                                                       92214, 92304, 140388,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 206776, 0, 3,
                                                                       204676, 138561, 204760,
                                                                       92304, 92394, 140514,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 206944, 0, 3,
                                                                       204760, 138624, 204844,
                                                                       92394, 92484, 140640,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 207112, 0, 3,
                                                                       204844, 138687, 204928,
                                                                       92484, 92574, 140766,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 207280, 0, 3,
                                                                       204928, 138750, 205012,
                                                                       92574, 92664, 140892,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 207448, 0, 3,
                                                                       205096, 139002, 205180,
                                                                       92844, 92934, 141270,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 207616, 0, 3,
                                                                       205180, 139065, 205264,
                                                                       92934, 93024, 141396,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 207784, 0, 3,
                                                                       205264, 139128, 205348,
                                                                       93024, 93114, 141522,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 207952, 0, 3,
                                                                       205348, 139191, 205432,
                                                                       93114, 93204, 141648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 208120, 0, 3,
                                                                       205432, 139254, 205516,
                                                                       93204, 93294, 141774,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 208288, 0, 3,
                                                                       205516, 139317, 205600,
                                                                       93294, 93384, 141900,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 208456, 0, 3,
                                                                       205600, 139380, 205684,
                                                                       93384, 93474, 142026,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 208624, 0, 3,
                                                                       205684, 139443, 205768,
                                                                       93474, 93564, 142152,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 208792, 0, 3,
                                                                       205768, 139506, 205852,
                                                                       93564, 93654, 142278,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 208960, 0, 3,
                                                                       205936, 139884, 206104,
                                                                       93834, 93984, 142824,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 209240, 0, 3,
                                                                       206104, 140010, 206272,
                                                                       93984, 94134, 143034,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 209520, 0, 3,
                                                                       206272, 140136, 206440,
                                                                       94134, 94284, 143244,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 209800, 0, 3,
                                                                       206440, 140262, 206608,
                                                                       94284, 94434, 143454,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 210080, 0, 3,
                                                                       206608, 140388, 206776,
                                                                       94434, 94584, 143664,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 210360, 0, 3,
                                                                       206776, 140514, 206944,
                                                                       94584, 94734, 143874,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 210640, 0, 3,
                                                                       206944, 140640, 207112,
                                                                       94734, 94884, 144084,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 210920, 0, 3,
                                                                       207112, 140766, 207280,
                                                                       94884, 95034, 144294,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 211200, 0, 3,
                                                                       207448, 141270, 207616,
                                                                       95334, 95484, 144924,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 211480, 0, 3,
                                                                       207616, 141396, 207784,
                                                                       95484, 95634, 145134,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 211760, 0, 3,
                                                                       207784, 141522, 207952,
                                                                       95634, 95784, 145344,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 212040, 0, 3,
                                                                       207952, 141648, 208120,
                                                                       95784, 95934, 145554,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 212320, 0, 3,
                                                                       208120, 141774, 208288,
                                                                       95934, 96084, 145764,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 212600, 0, 3,
                                                                       208288, 141900, 208456,
                                                                       96084, 96234, 145974,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 212880, 0, 3,
                                                                       208456, 142026, 208624,
                                                                       96234, 96384, 146184,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 213160, 0, 3,
                                                                       208624, 142152, 208792,
                                                                       96384, 96534, 146394,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 213440, 0, 3,
                                                                       208960, 142824, 209240,
                                                                       96834, 97059, 147234,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 213860, 0, 3,
                                                                       209240, 143034, 209520,
                                                                       97059, 97284, 147549,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 214280, 0, 3,
                                                                       209520, 143244, 209800,
                                                                       97284, 97509, 147864,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 214700, 0, 3,
                                                                       209800, 143454, 210080,
                                                                       97509, 97734, 148179,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 215120, 0, 3,
                                                                       210080, 143664, 210360,
                                                                       97734, 97959, 148494,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 215540, 0, 3,
                                                                       210360, 143874, 210640,
                                                                       97959, 98184, 148809,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 215960, 0, 3,
                                                                       210640, 144084, 210920,
                                                                       98184, 98409, 149124,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 216380, 0, 3,
                                                                       211200, 144924, 211480,
                                                                       98859, 99084, 150069,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 216800, 0, 3,
                                                                       211480, 145134, 211760,
                                                                       99084, 99309, 150384,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 217220, 0, 3,
                                                                       211760, 145344, 212040,
                                                                       99309, 99534, 150699,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 217640, 0, 3,
                                                                       212040, 145554, 212320,
                                                                       99534, 99759, 151014,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 218060, 0, 3,
                                                                       212320, 145764, 212600,
                                                                       99759, 99984, 151329,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 218480, 0, 3,
                                                                       212600, 145974, 212880,
                                                                       99984, 100209, 151644,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 218900, 0, 3,
                                                                       212880, 146184, 213160,
                                                                       100209, 100434, 151959,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 219320, 0, 3,
                                                                       213440, 147234, 213860,
                                                                       100884, 101199, 153156,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 219908, 0, 3,
                                                                       213860, 147549, 214280,
                                                                       101199, 101514, 153597,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 220496, 0, 3,
                                                                       214280, 147864, 214700,
                                                                       101514, 101829, 154038,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 221084, 0, 3,
                                                                       214700, 148179, 215120,
                                                                       101829, 102144, 154479,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 221672, 0, 3,
                                                                       215120, 148494, 215540,
                                                                       102144, 102459, 154920,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 222260, 0, 3,
                                                                       215540, 148809, 215960,
                                                                       102459, 102774, 155361,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 222848, 0, 3,
                                                                       216380, 150069, 216800,
                                                                       103404, 103719, 156684,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 223436, 0, 3,
                                                                       216800, 150384, 217220,
                                                                       103719, 104034, 157125,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 224024, 0, 3,
                                                                       217220, 150699, 217640,
                                                                       104034, 104349, 157566,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 224612, 0, 3,
                                                                       217640, 151014, 218060,
                                                                       104349, 104664, 158007,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 225200, 0, 3,
                                                                       218060, 151329, 218480,
                                                                       104664, 104979, 158448,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 225788, 0, 3,
                                                                       218480, 151644, 218900,
                                                                       104979, 105294, 158889,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 226376, 0, 3,
                                                                       219320, 153156, 219908,
                                                                       105924, 106344, 160506,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 227160, 0, 3,
                                                                       219908, 153597, 220496,
                                                                       106344, 106764, 161094,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 227944, 0, 3,
                                                                       220496, 154038, 221084,
                                                                       106764, 107184, 161682,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 228728, 0, 3,
                                                                       221084, 154479, 221672,
                                                                       107184, 107604, 162270,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 229512, 0, 3,
                                                                       221672, 154920, 222260,
                                                                       107604, 108024, 162858,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 230296, 0, 3,
                                                                       222848, 156684, 223436,
                                                                       108864, 109284, 164622,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 231080, 0, 3,
                                                                       223436, 157125, 224024,
                                                                       109284, 109704, 165210,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 231864, 0, 3,
                                                                       224024, 157566, 224612,
                                                                       109704, 110124, 165798,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 232648, 0, 3,
                                                                       224612, 158007, 225200,
                                                                       110124, 110544, 166386,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 233432, 0, 3,
                                                                       225200, 158448, 225788,
                                                                       110544, 110964, 166974,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 234216, 0, 3,
                                                                       226376, 160506, 227160,
                                                                       111804, 112344, 169074,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 235224, 0, 3,
                                                                       227160, 161094, 227944,
                                                                       112344, 112884, 169830,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 236232, 0, 3,
                                                                       227944, 161682, 228728,
                                                                       112884, 113424, 170586,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 237240, 0, 3,
                                                                       228728, 162270, 229512,
                                                                       113424, 113964, 171342,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 238248, 0, 3,
                                                                       230296, 164622, 231080,
                                                                       115044, 115584, 173610,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 239256, 0, 3,
                                                                       231080, 165210, 231864,
                                                                       115584, 116124, 174366,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 240264, 0, 3,
                                                                       231864, 165798, 232648,
                                                                       116124, 116664, 175122,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 241272, 0, 3,
                                                                       232648, 166386, 233432,
                                                                       116664, 117204, 175878,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 242280, 0, 3,
                                                                       234216, 169074, 235224,
                                                                       118284, 118959, 178524,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 243540, 0, 3,
                                                                       235224, 169830, 236232,
                                                                       118959, 119634, 179469,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 244800, 0, 3,
                                                                       236232, 170586, 237240,
                                                                       119634, 120309, 180414,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 246060, 0, 3,
                                                                       238248, 173610, 239256,
                                                                       121659, 122334, 183249,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 247320, 0, 3,
                                                                       239256, 174366, 240264,
                                                                       122334, 123009, 184194,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 248580, 0, 3,
                                                                       240264, 175122, 241272,
                                                                       123009, 123684, 185139,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 249840, 0, 3,
                                                                       242280, 178524, 243540,
                                                                       125034, 125859, 188394,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 251380, 0, 3,
                                                                       243540, 179469, 244800,
                                                                       125859, 126684, 189549,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 252920, 0, 3,
                                                                       246060, 183249, 247320,
                                                                       128334, 129159, 193014,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 254460, 0, 3,
                                                                       247320, 184194, 248580,
                                                                       129159, 129984, 194169,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 256000, 0, 3,
                                                                       249840, 188394, 251380,
                                                                       131634, 132624, 198096,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 257848, 0, 3,
                                                                       252920, 193014, 254460,
                                                                       134604, 135594, 202254,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259696, 3, 137574,
                                                                       137595, 203640, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259732, 3, 137595,
                                                                       137616, 203668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259768, 3, 137616,
                                                                       137637, 203696, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259804, 3, 137637,
                                                                       137658, 203724, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259840, 3, 137658,
                                                                       137679, 203752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259876, 3, 137679,
                                                                       137700, 203780, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259912, 3, 137700,
                                                                       137721, 203808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259948, 3, 137721,
                                                                       137742, 203836, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259984, 3, 137742,
                                                                       137763, 203864, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260020, 3, 137763,
                                                                       137784, 203892, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260056, 3, 137784,
                                                                       137805, 203920, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260092, 3, 137847,
                                                                       137868, 203948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260128, 3, 137868,
                                                                       137889, 203976, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260164, 3, 137889,
                                                                       137910, 204004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260200, 3, 137910,
                                                                       137931, 204032, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260236, 3, 137931,
                                                                       137952, 204060, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260272, 3, 137952,
                                                                       137973, 204088, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260308, 3, 137973,
                                                                       137994, 204116, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260344, 3, 137994,
                                                                       138015, 204144, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260380, 3, 138015,
                                                                       138036, 204172, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260416, 3, 138036,
                                                                       138057, 204200, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260452, 3, 138057,
                                                                       138078, 204228, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 260488, 0, 3,
                                                                       259696, 203640, 259732,
                                                                       138120, 138183, 204256,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 260596, 0, 3,
                                                                       259732, 203668, 259768,
                                                                       138183, 138246, 204340,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 260704, 0, 3,
                                                                       259768, 203696, 259804,
                                                                       138246, 138309, 204424,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 260812, 0, 3,
                                                                       259804, 203724, 259840,
                                                                       138309, 138372, 204508,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 260920, 0, 3,
                                                                       259840, 203752, 259876,
                                                                       138372, 138435, 204592,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 261028, 0, 3,
                                                                       259876, 203780, 259912,
                                                                       138435, 138498, 204676,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 261136, 0, 3,
                                                                       259912, 203808, 259948,
                                                                       138498, 138561, 204760,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 261244, 0, 3,
                                                                       259948, 203836, 259984,
                                                                       138561, 138624, 204844,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 261352, 0, 3,
                                                                       259984, 203864, 260020,
                                                                       138624, 138687, 204928,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 261460, 0, 3,
                                                                       260020, 203892, 260056,
                                                                       138687, 138750, 205012,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 261568, 0, 3,
                                                                       260092, 203948, 260128,
                                                                       138876, 138939, 205096,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 261676, 0, 3,
                                                                       260128, 203976, 260164,
                                                                       138939, 139002, 205180,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 261784, 0, 3,
                                                                       260164, 204004, 260200,
                                                                       139002, 139065, 205264,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 261892, 0, 3,
                                                                       260200, 204032, 260236,
                                                                       139065, 139128, 205348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 262000, 0, 3,
                                                                       260236, 204060, 260272,
                                                                       139128, 139191, 205432,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 262108, 0, 3,
                                                                       260272, 204088, 260308,
                                                                       139191, 139254, 205516,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 262216, 0, 3,
                                                                       260308, 204116, 260344,
                                                                       139254, 139317, 205600,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 262324, 0, 3,
                                                                       260344, 204144, 260380,
                                                                       139317, 139380, 205684,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 262432, 0, 3,
                                                                       260380, 204172, 260416,
                                                                       139380, 139443, 205768,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 262540, 0, 3,
                                                                       260416, 204200, 260452,
                                                                       139443, 139506, 205852,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 262648, 0, 3,
                                                                       260488, 204256, 260596,
                                                                       139632, 139758, 205936,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 262864, 0, 3,
                                                                       260596, 204340, 260704,
                                                                       139758, 139884, 206104,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 263080, 0, 3,
                                                                       260704, 204424, 260812,
                                                                       139884, 140010, 206272,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 263296, 0, 3,
                                                                       260812, 204508, 260920,
                                                                       140010, 140136, 206440,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 263512, 0, 3,
                                                                       260920, 204592, 261028,
                                                                       140136, 140262, 206608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 263728, 0, 3,
                                                                       261028, 204676, 261136,
                                                                       140262, 140388, 206776,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 263944, 0, 3,
                                                                       261136, 204760, 261244,
                                                                       140388, 140514, 206944,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 264160, 0, 3,
                                                                       261244, 204844, 261352,
                                                                       140514, 140640, 207112,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 264376, 0, 3,
                                                                       261352, 204928, 261460,
                                                                       140640, 140766, 207280,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 264592, 0, 3,
                                                                       261568, 205096, 261676,
                                                                       141018, 141144, 207448,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 264808, 0, 3,
                                                                       261676, 205180, 261784,
                                                                       141144, 141270, 207616,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 265024, 0, 3,
                                                                       261784, 205264, 261892,
                                                                       141270, 141396, 207784,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 265240, 0, 3,
                                                                       261892, 205348, 262000,
                                                                       141396, 141522, 207952,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 265456, 0, 3,
                                                                       262000, 205432, 262108,
                                                                       141522, 141648, 208120,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 265672, 0, 3,
                                                                       262108, 205516, 262216,
                                                                       141648, 141774, 208288,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 265888, 0, 3,
                                                                       262216, 205600, 262324,
                                                                       141774, 141900, 208456,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 266104, 0, 3,
                                                                       262324, 205684, 262432,
                                                                       141900, 142026, 208624,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 266320, 0, 3,
                                                                       262432, 205768, 262540,
                                                                       142026, 142152, 208792,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 266536, 0, 3,
                                                                       262648, 205936, 262864,
                                                                       142404, 142614, 208960,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 266896, 0, 3,
                                                                       262864, 206104, 263080,
                                                                       142614, 142824, 209240,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 267256, 0, 3,
                                                                       263080, 206272, 263296,
                                                                       142824, 143034, 209520,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 267616, 0, 3,
                                                                       263296, 206440, 263512,
                                                                       143034, 143244, 209800,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 267976, 0, 3,
                                                                       263512, 206608, 263728,
                                                                       143244, 143454, 210080,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 268336, 0, 3,
                                                                       263728, 206776, 263944,
                                                                       143454, 143664, 210360,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 268696, 0, 3,
                                                                       263944, 206944, 264160,
                                                                       143664, 143874, 210640,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 269056, 0, 3,
                                                                       264160, 207112, 264376,
                                                                       143874, 144084, 210920,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 269416, 0, 3,
                                                                       264592, 207448, 264808,
                                                                       144504, 144714, 211200,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 269776, 0, 3,
                                                                       264808, 207616, 265024,
                                                                       144714, 144924, 211480,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 270136, 0, 3,
                                                                       265024, 207784, 265240,
                                                                       144924, 145134, 211760,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 270496, 0, 3,
                                                                       265240, 207952, 265456,
                                                                       145134, 145344, 212040,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 270856, 0, 3,
                                                                       265456, 208120, 265672,
                                                                       145344, 145554, 212320,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 271216, 0, 3,
                                                                       265672, 208288, 265888,
                                                                       145554, 145764, 212600,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 271576, 0, 3,
                                                                       265888, 208456, 266104,
                                                                       145764, 145974, 212880,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 271936, 0, 3,
                                                                       266104, 208624, 266320,
                                                                       145974, 146184, 213160,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 272296, 0, 3,
                                                                       266536, 208960, 266896,
                                                                       146604, 146919, 213440,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 272836, 0, 3,
                                                                       266896, 209240, 267256,
                                                                       146919, 147234, 213860,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 273376, 0, 3,
                                                                       267256, 209520, 267616,
                                                                       147234, 147549, 214280,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 273916, 0, 3,
                                                                       267616, 209800, 267976,
                                                                       147549, 147864, 214700,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 274456, 0, 3,
                                                                       267976, 210080, 268336,
                                                                       147864, 148179, 215120,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 274996, 0, 3,
                                                                       268336, 210360, 268696,
                                                                       148179, 148494, 215540,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 275536, 0, 3,
                                                                       268696, 210640, 269056,
                                                                       148494, 148809, 215960,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 276076, 0, 3,
                                                                       269416, 211200, 269776,
                                                                       149439, 149754, 216380,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 276616, 0, 3,
                                                                       269776, 211480, 270136,
                                                                       149754, 150069, 216800,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 277156, 0, 3,
                                                                       270136, 211760, 270496,
                                                                       150069, 150384, 217220,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 277696, 0, 3,
                                                                       270496, 212040, 270856,
                                                                       150384, 150699, 217640,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 278236, 0, 3,
                                                                       270856, 212320, 271216,
                                                                       150699, 151014, 218060,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 278776, 0, 3,
                                                                       271216, 212600, 271576,
                                                                       151014, 151329, 218480,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 279316, 0, 3,
                                                                       271576, 212880, 271936,
                                                                       151329, 151644, 218900,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 279856, 0, 3,
                                                                       272296, 213440, 272836,
                                                                       152274, 152715, 219320,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 280612, 0, 3,
                                                                       272836, 213860, 273376,
                                                                       152715, 153156, 219908,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 281368, 0, 3,
                                                                       273376, 214280, 273916,
                                                                       153156, 153597, 220496,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 282124, 0, 3,
                                                                       273916, 214700, 274456,
                                                                       153597, 154038, 221084,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 282880, 0, 3,
                                                                       274456, 215120, 274996,
                                                                       154038, 154479, 221672,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 283636, 0, 3,
                                                                       274996, 215540, 275536,
                                                                       154479, 154920, 222260,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 284392, 0, 3,
                                                                       276076, 216380, 276616,
                                                                       155802, 156243, 222848,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 285148, 0, 3,
                                                                       276616, 216800, 277156,
                                                                       156243, 156684, 223436,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 285904, 0, 3,
                                                                       277156, 217220, 277696,
                                                                       156684, 157125, 224024,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 286660, 0, 3,
                                                                       277696, 217640, 278236,
                                                                       157125, 157566, 224612,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 287416, 0, 3,
                                                                       278236, 218060, 278776,
                                                                       157566, 158007, 225200,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 288172, 0, 3,
                                                                       278776, 218480, 279316,
                                                                       158007, 158448, 225788,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 288928, 0, 3,
                                                                       279856, 219320, 280612,
                                                                       159330, 159918, 226376,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 289936, 0, 3,
                                                                       280612, 219908, 281368,
                                                                       159918, 160506, 227160,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 290944, 0, 3,
                                                                       281368, 220496, 282124,
                                                                       160506, 161094, 227944,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 291952, 0, 3,
                                                                       282124, 221084, 282880,
                                                                       161094, 161682, 228728,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 292960, 0, 3,
                                                                       282880, 221672, 283636,
                                                                       161682, 162270, 229512,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 293968, 0, 3,
                                                                       284392, 222848, 285148,
                                                                       163446, 164034, 230296,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 294976, 0, 3,
                                                                       285148, 223436, 285904,
                                                                       164034, 164622, 231080,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 295984, 0, 3,
                                                                       285904, 224024, 286660,
                                                                       164622, 165210, 231864,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 296992, 0, 3,
                                                                       286660, 224612, 287416,
                                                                       165210, 165798, 232648,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 298000, 0, 3,
                                                                       287416, 225200, 288172,
                                                                       165798, 166386, 233432,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 299008, 0, 3,
                                                                       288928, 226376, 289936,
                                                                       167562, 168318, 234216,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 300304, 0, 3,
                                                                       289936, 227160, 290944,
                                                                       168318, 169074, 235224,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 301600, 0, 3,
                                                                       290944, 227944, 291952,
                                                                       169074, 169830, 236232,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 302896, 0, 3,
                                                                       291952, 228728, 292960,
                                                                       169830, 170586, 237240,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 304192, 0, 3,
                                                                       293968, 230296, 294976,
                                                                       172098, 172854, 238248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 305488, 0, 3,
                                                                       294976, 231080, 295984,
                                                                       172854, 173610, 239256,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 306784, 0, 3,
                                                                       295984, 231864, 296992,
                                                                       173610, 174366, 240264,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 308080, 0, 3,
                                                                       296992, 232648, 298000,
                                                                       174366, 175122, 241272,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 309376, 0, 3,
                                                                       299008, 234216, 300304,
                                                                       176634, 177579, 242280,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 310996, 0, 3,
                                                                       300304, 235224, 301600,
                                                                       177579, 178524, 243540,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 312616, 0, 3,
                                                                       301600, 236232, 302896,
                                                                       178524, 179469, 244800,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 314236, 0, 3,
                                                                       304192, 238248, 305488,
                                                                       181359, 182304, 246060,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 315856, 0, 3,
                                                                       305488, 239256, 306784,
                                                                       182304, 183249, 247320,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 317476, 0, 3,
                                                                       306784, 240264, 308080,
                                                                       183249, 184194, 248580,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 319096, 0, 3,
                                                                       309376, 242280, 310996,
                                                                       186084, 187239, 249840,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 321076, 0, 3,
                                                                       310996, 243540, 312616,
                                                                       187239, 188394, 251380,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 323056, 0, 3,
                                                                       314236, 246060, 315856,
                                                                       190704, 191859, 252920,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 325036, 0, 3,
                                                                       315856, 247320, 317476,
                                                                       191859, 193014, 254460,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 327016, 0, 3,
                                                                       319096, 249840, 321076,
                                                                       195324, 196710, 256000,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 329392, 0, 3,
                                                                       323056, 252920, 325036,
                                                                       199482, 200868, 257848,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 331768, 288928, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 333196, 293968, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 334624, 299008, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 336460, 304192, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 338296, 309376, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 340591, 314236, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 342886, 319096, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 345691, 323056, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 348496, 327016, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 351862, 329392, 2376, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 332776, 331768, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 334204, 333196, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 335920, 334624, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 337756, 336460, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 339916, 338296, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 342211, 340591, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 344866, 342886, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 347671, 345691, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 350872, 348496, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 354238, 351862, 66, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 355228, 332776, 335920, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 356488, 334204, 337756, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 357748, 335920, 339916, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 359368, 337756, 342211, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 360988, 339916, 344866, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 363013, 342211, 347671, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 365038, 344866, 350872, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 367513, 347671, 354238, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 369988, 355228, 357748, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 372508, 356488, 359368, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 375028, 357748, 360988, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 378268, 359368, 363013, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 381508, 360988, 365038, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 385558, 363013, 367513, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 389608, 369988, 375028, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 393808, 372508, 378268, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 398008, 375028, 381508, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 403408, 378268, 385558, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 408808, 389608, 398008, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 415108, 393808, 403408, 15,
                                             nmax);

        simdtrf::transform_g_inner(buffer, 421408, 415108, 28, 15, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 421408, 135, nmax);

        simdtrf::transform_g_inner(buffer, 421408, 408808, 28, 15, nmax);

        simdtrf::transform_i_outer(values + 1755 * nvalues + n * npairs, nvalues, buffer, 421408,
                                   135, nmax);
    }

    for (size_t m = 0; m < 3510; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
