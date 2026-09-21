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


#include "SimdThreeCenterElectronRepulsionRsRecIIH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
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
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_iih_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_iih_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 446622, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3718 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 446622, 247528, 24535, dimensions);

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

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5442, 0, 3, 3638,
                                                                       3693, 4518, 4584, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5520, 0, 3, 3693,
                                                                       3748, 4584, 4650, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5598, 0, 3, 3748,
                                                                       3803, 4650, 4716, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5676, 0, 3, 3803,
                                                                       3858, 4716, 4782, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5754, 0, 3, 3858,
                                                                       3913, 4782, 4848, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5832, 0, 3, 3913,
                                                                       3968, 4848, 4914, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5910, 0, 3, 4078,
                                                                       4133, 4980, 5046, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5988, 0, 3, 4133,
                                                                       4188, 5046, 5112, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6066, 0, 3, 4188,
                                                                       4243, 5112, 5178, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6144, 0, 3, 4243,
                                                                       4298, 5178, 5244, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6222, 0, 3, 4298,
                                                                       4353, 5244, 5310, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6300, 0, 3, 4353,
                                                                       4408, 5310, 5376, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 6378, 0, 3, 4518,
                                                                       4584, 5442, 5520, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 6469, 0, 3, 4584,
                                                                       4650, 5520, 5598, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 6560, 0, 3, 4650,
                                                                       4716, 5598, 5676, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 6651, 0, 3, 4716,
                                                                       4782, 5676, 5754, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 6742, 0, 3, 4782,
                                                                       4848, 5754, 5832, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 6833, 0, 3, 4980,
                                                                       5046, 5910, 5988, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 6924, 0, 3, 5046,
                                                                       5112, 5988, 6066, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 7015, 0, 3, 5112,
                                                                       5178, 6066, 6144, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 7106, 0, 3, 5178,
                                                                       5244, 6144, 6222, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 7197, 0, 3, 5244,
                                                                       5310, 6222, 6300, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7288, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7291, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7294, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7297, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7300, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7303, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7306, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7309, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7312, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7315, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7318, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7321, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7324, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7327, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7330, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7333, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7336, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7339, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7342, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7345, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7348, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7351, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7354, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7357, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7360, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7363, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7366, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7369, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7372, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7375, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7378, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7381, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7384, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7387, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7390, 3, 7, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7399, 3, 8, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7408, 3, 9, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7417, 3, 10, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7426, 3, 11, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7435, 3, 12, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7444, 3, 13, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7453, 3, 14, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7462, 3, 15, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7471, 3, 16, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7480, 3, 17, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7489, 3, 18, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7498, 3, 19, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7507, 3, 20, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7516, 3, 21, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7525, 3, 22, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7534, 3, 25, 90,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7543, 3, 26, 93,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7552, 3, 27, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7561, 3, 28, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7570, 3, 29, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7579, 3, 30, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7588, 3, 31, 108,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7597, 3, 32, 111,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7606, 3, 33, 114,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7615, 3, 34, 117,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7624, 3, 35, 120,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7633, 3, 36, 123,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7642, 3, 37, 126,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7651, 3, 38, 129,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7660, 3, 39, 132,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7669, 3, 40, 135,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7678, 3, 42, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7696, 3, 45, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7714, 3, 48, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7732, 3, 51, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7750, 3, 54, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7768, 3, 57, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7786, 3, 60, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7804, 3, 63, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7822, 3, 66, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7840, 3, 69, 192,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7858, 3, 72, 198,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7876, 3, 75, 204,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7894, 3, 78, 210,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7912, 3, 81, 216,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7930, 3, 84, 222,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7948, 3, 90, 228,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7966, 3, 93, 234,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7984, 3, 96, 240,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8002, 3, 99, 246,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8020, 3, 102, 252,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8038, 3, 105, 258,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8056, 3, 108, 264,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8074, 3, 111, 270,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8092, 3, 114, 276,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8110, 3, 117, 282,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8128, 3, 120, 288,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8146, 3, 123, 294,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8164, 3, 126, 300,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8182, 3, 129, 306,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8200, 3, 132, 312,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8218, 3, 138, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8248, 3, 144, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8278, 3, 150, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8308, 3, 156, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8338, 3, 162, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8368, 3, 168, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8398, 3, 174, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8428, 3, 180, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8458, 3, 186, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8488, 3, 192, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8518, 3, 198, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8548, 3, 204, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8578, 3, 210, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8608, 3, 216, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8638, 3, 228, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8668, 3, 234, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8698, 3, 240, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8728, 3, 246, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8758, 3, 252, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8788, 3, 258, 508,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8818, 3, 264, 518,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8848, 3, 270, 528,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8878, 3, 276, 538,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8908, 3, 282, 548,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8938, 3, 288, 558,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8968, 3, 294, 568,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8998, 3, 300, 578,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 9028, 3, 306, 588,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9058, 3, 318, 598,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9103, 3, 328, 613,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9148, 3, 338, 628,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9193, 3, 348, 643,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9238, 3, 358, 658,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9283, 3, 368, 673,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9328, 3, 378, 688,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9373, 3, 388, 703,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9418, 3, 398, 718,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9463, 3, 408, 733,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9508, 3, 418, 748,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9553, 3, 428, 763,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9598, 3, 438, 778,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9643, 3, 458, 793,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9688, 3, 468, 808,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9733, 3, 478, 823,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9778, 3, 488, 838,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9823, 3, 498, 853,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9868, 3, 508, 868,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9913, 3, 518, 883,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9958, 3, 528, 898,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10003, 3, 538,
                                                                       913, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10048, 3, 548,
                                                                       928, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10093, 3, 558,
                                                                       943, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10138, 3, 568,
                                                                       958, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10183, 3, 578,
                                                                       973, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10228, 3, 598,
                                                                       988, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10291, 3, 613,
                                                                       1009, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10354, 3, 628,
                                                                       1030, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10417, 3, 643,
                                                                       1051, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10480, 3, 658,
                                                                       1072, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10543, 3, 673,
                                                                       1093, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10606, 3, 688,
                                                                       1114, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10669, 3, 703,
                                                                       1135, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10732, 3, 718,
                                                                       1156, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10795, 3, 733,
                                                                       1177, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10858, 3, 748,
                                                                       1198, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10921, 3, 763,
                                                                       1219, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10984, 3, 793,
                                                                       1240, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11047, 3, 808,
                                                                       1261, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11110, 3, 823,
                                                                       1282, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11173, 3, 838,
                                                                       1303, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11236, 3, 853,
                                                                       1324, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11299, 3, 868,
                                                                       1345, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11362, 3, 883,
                                                                       1366, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11425, 3, 898,
                                                                       1387, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11488, 3, 913,
                                                                       1408, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11551, 3, 928,
                                                                       1429, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11614, 3, 943,
                                                                       1450, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11677, 3, 958,
                                                                       1471, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11740, 3, 988,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11824, 3, 1009,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11908, 3, 1030,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11992, 3, 1051,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12076, 3, 1072,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12160, 3, 1093,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12244, 3, 1114,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12328, 3, 1135,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12412, 3, 1156,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12496, 3, 1177,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12580, 3, 1198,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12664, 3, 1240,
                                                                       1800, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12748, 3, 1261,
                                                                       1828, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12832, 3, 1282,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12916, 3, 1303,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13000, 3, 1324,
                                                                       1912, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13084, 3, 1345,
                                                                       1940, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13168, 3, 1366,
                                                                       1968, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13252, 3, 1387,
                                                                       1996, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13336, 3, 1408,
                                                                       2024, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13420, 3, 1429,
                                                                       2052, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13504, 3, 1450,
                                                                       2080, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13588, 3, 1492,
                                                                       2108, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13696, 3, 1520,
                                                                       2144, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13804, 3, 1548,
                                                                       2180, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13912, 3, 1576,
                                                                       2216, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14020, 3, 1604,
                                                                       2252, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14128, 3, 1632,
                                                                       2288, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14236, 3, 1660,
                                                                       2324, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14344, 3, 1688,
                                                                       2360, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14452, 3, 1716,
                                                                       2396, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14560, 3, 1744,
                                                                       2432, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14668, 3, 1800,
                                                                       2468, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14776, 3, 1828,
                                                                       2504, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14884, 3, 1856,
                                                                       2540, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14992, 3, 1884,
                                                                       2576, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15100, 3, 1912,
                                                                       2612, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15208, 3, 1940,
                                                                       2648, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15316, 3, 1968,
                                                                       2684, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15424, 3, 1996,
                                                                       2720, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15532, 3, 2024,
                                                                       2756, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15640, 3, 2052,
                                                                       2792, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15748, 3, 2108,
                                                                       2828, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15883, 3, 2144,
                                                                       2873, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16018, 3, 2180,
                                                                       2918, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16153, 3, 2216,
                                                                       2963, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16288, 3, 2252,
                                                                       3008, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16423, 3, 2288,
                                                                       3053, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16558, 3, 2324,
                                                                       3098, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16693, 3, 2360,
                                                                       3143, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16828, 3, 2396,
                                                                       3188, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16963, 3, 2468,
                                                                       3233, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17098, 3, 2504,
                                                                       3278, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17233, 3, 2540,
                                                                       3323, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17368, 3, 2576,
                                                                       3368, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17503, 3, 2612,
                                                                       3413, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17638, 3, 2648,
                                                                       3458, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17773, 3, 2684,
                                                                       3503, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17908, 3, 2720,
                                                                       3548, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18043, 3, 2756,
                                                                       3593, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18178, 3, 2828,
                                                                       3638, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18343, 3, 2873,
                                                                       3693, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18508, 3, 2918,
                                                                       3748, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18673, 3, 2963,
                                                                       3803, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18838, 3, 3008,
                                                                       3858, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19003, 3, 3053,
                                                                       3913, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19168, 3, 3098,
                                                                       3968, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19333, 3, 3143,
                                                                       4023, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19498, 3, 3233,
                                                                       4078, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19663, 3, 3278,
                                                                       4133, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19828, 3, 3323,
                                                                       4188, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19993, 3, 3368,
                                                                       4243, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20158, 3, 3413,
                                                                       4298, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20323, 3, 3458,
                                                                       4353, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20488, 3, 3503,
                                                                       4408, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20653, 3, 3548,
                                                                       4463, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20818, 3, 3638,
                                                                       4518, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21016, 3, 3693,
                                                                       4584, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21214, 3, 3748,
                                                                       4650, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21412, 3, 3803,
                                                                       4716, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21610, 3, 3858,
                                                                       4782, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21808, 3, 3913,
                                                                       4848, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22006, 3, 3968,
                                                                       4914, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22204, 3, 4078,
                                                                       4980, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22402, 3, 4133,
                                                                       5046, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22600, 3, 4188,
                                                                       5112, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22798, 3, 4243,
                                                                       5178, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22996, 3, 4298,
                                                                       5244, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23194, 3, 4353,
                                                                       5310, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23392, 3, 4408,
                                                                       5376, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 23590, 3, 4518,
                                                                       5442, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 23824, 3, 4584,
                                                                       5520, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 24058, 3, 4650,
                                                                       5598, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 24292, 3, 4716,
                                                                       5676, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 24526, 3, 4782,
                                                                       5754, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 24760, 3, 4848,
                                                                       5832, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 24994, 3, 4980,
                                                                       5910, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 25228, 3, 5046,
                                                                       5988, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 25462, 3, 5112,
                                                                       6066, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 25696, 3, 5178,
                                                                       6144, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 25930, 3, 5244,
                                                                       6222, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 26164, 3, 5310,
                                                                       6300, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 26398, 3, 5442,
                                                                       6378, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 26671, 3, 5520,
                                                                       6469, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 26944, 3, 5598,
                                                                       6560, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 27217, 3, 5676,
                                                                       6651, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 27490, 3, 5754,
                                                                       6742, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 27763, 3, 5910,
                                                                       6833, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 28036, 3, 5988,
                                                                       6924, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 28309, 3, 6066,
                                                                       7015, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 28582, 3, 6144,
                                                                       7106, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 28855, 3, 6222,
                                                                       7197, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29128, 3, 7, 8,
                                                                       7294, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29134, 3, 8, 9,
                                                                       7297, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29140, 3, 9, 10,
                                                                       7300, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29146, 3, 10, 11,
                                                                       7303, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29152, 3, 11, 12,
                                                                       7306, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29158, 3, 12, 13,
                                                                       7309, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29164, 3, 13, 14,
                                                                       7312, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29170, 3, 14, 15,
                                                                       7315, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29176, 3, 15, 16,
                                                                       7318, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29182, 3, 16, 17,
                                                                       7321, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29188, 3, 17, 18,
                                                                       7324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29194, 3, 18, 19,
                                                                       7327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29200, 3, 19, 20,
                                                                       7330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29206, 3, 20, 21,
                                                                       7333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29212, 3, 21, 22,
                                                                       7336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29218, 3, 25, 26,
                                                                       7345, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29224, 3, 26, 27,
                                                                       7348, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29230, 3, 27, 28,
                                                                       7351, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29236, 3, 28, 29,
                                                                       7354, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29242, 3, 29, 30,
                                                                       7357, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29248, 3, 30, 31,
                                                                       7360, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29254, 3, 31, 32,
                                                                       7363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29260, 3, 32, 33,
                                                                       7366, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29266, 3, 33, 34,
                                                                       7369, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29272, 3, 34, 35,
                                                                       7372, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29278, 3, 35, 36,
                                                                       7375, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29284, 3, 36, 37,
                                                                       7378, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29290, 3, 37, 38,
                                                                       7381, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29296, 3, 38, 39,
                                                                       7384, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 29302, 3, 39, 40,
                                                                       7387, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29308, 0, 3,
                                                                       29128, 7294, 29134, 7408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29326, 0, 3,
                                                                       29134, 7297, 29140, 7417,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29344, 0, 3,
                                                                       29140, 7300, 29146, 7426,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29362, 0, 3,
                                                                       29146, 7303, 29152, 7435,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29380, 0, 3,
                                                                       29152, 7306, 29158, 7444,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29398, 0, 3,
                                                                       29158, 7309, 29164, 7453,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29416, 0, 3,
                                                                       29164, 7312, 29170, 7462,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29434, 0, 3,
                                                                       29170, 7315, 29176, 7471,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29452, 0, 3,
                                                                       29176, 7318, 29182, 7480,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29470, 0, 3,
                                                                       29182, 7321, 29188, 7489,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29488, 0, 3,
                                                                       29188, 7324, 29194, 7498,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29506, 0, 3,
                                                                       29194, 7327, 29200, 7507,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29524, 0, 3,
                                                                       29200, 7330, 29206, 7516,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29542, 0, 3,
                                                                       29206, 7333, 29212, 7525,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29560, 0, 3,
                                                                       29218, 7345, 29224, 7552,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29578, 0, 3,
                                                                       29224, 7348, 29230, 7561,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29596, 0, 3,
                                                                       29230, 7351, 29236, 7570,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29614, 0, 3,
                                                                       29236, 7354, 29242, 7579,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29632, 0, 3,
                                                                       29242, 7357, 29248, 7588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29650, 0, 3,
                                                                       29248, 7360, 29254, 7597,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29668, 0, 3,
                                                                       29254, 7363, 29260, 7606,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29686, 0, 3,
                                                                       29260, 7366, 29266, 7615,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29704, 0, 3,
                                                                       29266, 7369, 29272, 7624,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29722, 0, 3,
                                                                       29272, 7372, 29278, 7633,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29740, 0, 3,
                                                                       29278, 7375, 29284, 7642,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29758, 0, 3,
                                                                       29284, 7378, 29290, 7651,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29776, 0, 3,
                                                                       29290, 7381, 29296, 7660,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29794, 0, 3,
                                                                       29296, 7384, 29302, 7669,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29812, 0, 3,
                                                                       29308, 7408, 29326, 138,
                                                                       144, 7714, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29848, 0, 3,
                                                                       29326, 7417, 29344, 144,
                                                                       150, 7732, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29884, 0, 3,
                                                                       29344, 7426, 29362, 150,
                                                                       156, 7750, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29920, 0, 3,
                                                                       29362, 7435, 29380, 156,
                                                                       162, 7768, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29956, 0, 3,
                                                                       29380, 7444, 29398, 162,
                                                                       168, 7786, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29992, 0, 3,
                                                                       29398, 7453, 29416, 168,
                                                                       174, 7804, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30028, 0, 3,
                                                                       29416, 7462, 29434, 174,
                                                                       180, 7822, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30064, 0, 3,
                                                                       29434, 7471, 29452, 180,
                                                                       186, 7840, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30100, 0, 3,
                                                                       29452, 7480, 29470, 186,
                                                                       192, 7858, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30136, 0, 3,
                                                                       29470, 7489, 29488, 192,
                                                                       198, 7876, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30172, 0, 3,
                                                                       29488, 7498, 29506, 198,
                                                                       204, 7894, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30208, 0, 3,
                                                                       29506, 7507, 29524, 204,
                                                                       210, 7912, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30244, 0, 3,
                                                                       29524, 7516, 29542, 210,
                                                                       216, 7930, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30280, 0, 3,
                                                                       29560, 7552, 29578, 228,
                                                                       234, 7984, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30316, 0, 3,
                                                                       29578, 7561, 29596, 234,
                                                                       240, 8002, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30352, 0, 3,
                                                                       29596, 7570, 29614, 240,
                                                                       246, 8020, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30388, 0, 3,
                                                                       29614, 7579, 29632, 246,
                                                                       252, 8038, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30424, 0, 3,
                                                                       29632, 7588, 29650, 252,
                                                                       258, 8056, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30460, 0, 3,
                                                                       29650, 7597, 29668, 258,
                                                                       264, 8074, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30496, 0, 3,
                                                                       29668, 7606, 29686, 264,
                                                                       270, 8092, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30532, 0, 3,
                                                                       29686, 7615, 29704, 270,
                                                                       276, 8110, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30568, 0, 3,
                                                                       29704, 7624, 29722, 276,
                                                                       282, 8128, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30604, 0, 3,
                                                                       29722, 7633, 29740, 282,
                                                                       288, 8146, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30640, 0, 3,
                                                                       29740, 7642, 29758, 288,
                                                                       294, 8164, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30676, 0, 3,
                                                                       29758, 7651, 29776, 294,
                                                                       300, 8182, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30712, 0, 3,
                                                                       29776, 7660, 29794, 300,
                                                                       306, 8200, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30748, 0, 3,
                                                                       29812, 7714, 29848, 318,
                                                                       328, 8278, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30808, 0, 3,
                                                                       29848, 7732, 29884, 328,
                                                                       338, 8308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30868, 0, 3,
                                                                       29884, 7750, 29920, 338,
                                                                       348, 8338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30928, 0, 3,
                                                                       29920, 7768, 29956, 348,
                                                                       358, 8368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30988, 0, 3,
                                                                       29956, 7786, 29992, 358,
                                                                       368, 8398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31048, 0, 3,
                                                                       29992, 7804, 30028, 368,
                                                                       378, 8428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31108, 0, 3,
                                                                       30028, 7822, 30064, 378,
                                                                       388, 8458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31168, 0, 3,
                                                                       30064, 7840, 30100, 388,
                                                                       398, 8488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31228, 0, 3,
                                                                       30100, 7858, 30136, 398,
                                                                       408, 8518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31288, 0, 3,
                                                                       30136, 7876, 30172, 408,
                                                                       418, 8548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31348, 0, 3,
                                                                       30172, 7894, 30208, 418,
                                                                       428, 8578, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31408, 0, 3,
                                                                       30208, 7912, 30244, 428,
                                                                       438, 8608, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31468, 0, 3,
                                                                       30280, 7984, 30316, 458,
                                                                       468, 8698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31528, 0, 3,
                                                                       30316, 8002, 30352, 468,
                                                                       478, 8728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31588, 0, 3,
                                                                       30352, 8020, 30388, 478,
                                                                       488, 8758, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31648, 0, 3,
                                                                       30388, 8038, 30424, 488,
                                                                       498, 8788, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31708, 0, 3,
                                                                       30424, 8056, 30460, 498,
                                                                       508, 8818, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31768, 0, 3,
                                                                       30460, 8074, 30496, 508,
                                                                       518, 8848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31828, 0, 3,
                                                                       30496, 8092, 30532, 518,
                                                                       528, 8878, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31888, 0, 3,
                                                                       30532, 8110, 30568, 528,
                                                                       538, 8908, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31948, 0, 3,
                                                                       30568, 8128, 30604, 538,
                                                                       548, 8938, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 32008, 0, 3,
                                                                       30604, 8146, 30640, 548,
                                                                       558, 8968, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 32068, 0, 3,
                                                                       30640, 8164, 30676, 558,
                                                                       568, 8998, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 32128, 0, 3,
                                                                       30676, 8182, 30712, 568,
                                                                       578, 9028, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32188, 0, 3,
                                                                       30748, 8278, 30808, 598,
                                                                       613, 9148, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32278, 0, 3,
                                                                       30808, 8308, 30868, 613,
                                                                       628, 9193, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32368, 0, 3,
                                                                       30868, 8338, 30928, 628,
                                                                       643, 9238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32458, 0, 3,
                                                                       30928, 8368, 30988, 643,
                                                                       658, 9283, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32548, 0, 3,
                                                                       30988, 8398, 31048, 658,
                                                                       673, 9328, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32638, 0, 3,
                                                                       31048, 8428, 31108, 673,
                                                                       688, 9373, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32728, 0, 3,
                                                                       31108, 8458, 31168, 688,
                                                                       703, 9418, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32818, 0, 3,
                                                                       31168, 8488, 31228, 703,
                                                                       718, 9463, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32908, 0, 3,
                                                                       31228, 8518, 31288, 718,
                                                                       733, 9508, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32998, 0, 3,
                                                                       31288, 8548, 31348, 733,
                                                                       748, 9553, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33088, 0, 3,
                                                                       31348, 8578, 31408, 748,
                                                                       763, 9598, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33178, 0, 3,
                                                                       31468, 8698, 31528, 793,
                                                                       808, 9733, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33268, 0, 3,
                                                                       31528, 8728, 31588, 808,
                                                                       823, 9778, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33358, 0, 3,
                                                                       31588, 8758, 31648, 823,
                                                                       838, 9823, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33448, 0, 3,
                                                                       31648, 8788, 31708, 838,
                                                                       853, 9868, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33538, 0, 3,
                                                                       31708, 8818, 31768, 853,
                                                                       868, 9913, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33628, 0, 3,
                                                                       31768, 8848, 31828, 868,
                                                                       883, 9958, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33718, 0, 3,
                                                                       31828, 8878, 31888, 883,
                                                                       898, 10003, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33808, 0, 3,
                                                                       31888, 8908, 31948, 898,
                                                                       913, 10048, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33898, 0, 3,
                                                                       31948, 8938, 32008, 913,
                                                                       928, 10093, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33988, 0, 3,
                                                                       32008, 8968, 32068, 928,
                                                                       943, 10138, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 34078, 0, 3,
                                                                       32068, 8998, 32128, 943,
                                                                       958, 10183, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34168, 0, 3,
                                                                       32188, 9148, 32278, 988,
                                                                       1009, 10354, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34294, 0, 3,
                                                                       32278, 9193, 32368, 1009,
                                                                       1030, 10417, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34420, 0, 3,
                                                                       32368, 9238, 32458, 1030,
                                                                       1051, 10480, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34546, 0, 3,
                                                                       32458, 9283, 32548, 1051,
                                                                       1072, 10543, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34672, 0, 3,
                                                                       32548, 9328, 32638, 1072,
                                                                       1093, 10606, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34798, 0, 3,
                                                                       32638, 9373, 32728, 1093,
                                                                       1114, 10669, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34924, 0, 3,
                                                                       32728, 9418, 32818, 1114,
                                                                       1135, 10732, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35050, 0, 3,
                                                                       32818, 9463, 32908, 1135,
                                                                       1156, 10795, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35176, 0, 3,
                                                                       32908, 9508, 32998, 1156,
                                                                       1177, 10858, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35302, 0, 3,
                                                                       32998, 9553, 33088, 1177,
                                                                       1198, 10921, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35428, 0, 3,
                                                                       33178, 9733, 33268, 1240,
                                                                       1261, 11110, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35554, 0, 3,
                                                                       33268, 9778, 33358, 1261,
                                                                       1282, 11173, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35680, 0, 3,
                                                                       33358, 9823, 33448, 1282,
                                                                       1303, 11236, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35806, 0, 3,
                                                                       33448, 9868, 33538, 1303,
                                                                       1324, 11299, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35932, 0, 3,
                                                                       33538, 9913, 33628, 1324,
                                                                       1345, 11362, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36058, 0, 3,
                                                                       33628, 9958, 33718, 1345,
                                                                       1366, 11425, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36184, 0, 3,
                                                                       33718, 10003, 33808, 1366,
                                                                       1387, 11488, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36310, 0, 3,
                                                                       33808, 10048, 33898, 1387,
                                                                       1408, 11551, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36436, 0, 3,
                                                                       33898, 10093, 33988, 1408,
                                                                       1429, 11614, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36562, 0, 3,
                                                                       33988, 10138, 34078, 1429,
                                                                       1450, 11677, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 36688, 0, 3,
                                                                       34168, 10354, 34294, 1492,
                                                                       1520, 11908, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 36856, 0, 3,
                                                                       34294, 10417, 34420, 1520,
                                                                       1548, 11992, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37024, 0, 3,
                                                                       34420, 10480, 34546, 1548,
                                                                       1576, 12076, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37192, 0, 3,
                                                                       34546, 10543, 34672, 1576,
                                                                       1604, 12160, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37360, 0, 3,
                                                                       34672, 10606, 34798, 1604,
                                                                       1632, 12244, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37528, 0, 3,
                                                                       34798, 10669, 34924, 1632,
                                                                       1660, 12328, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37696, 0, 3,
                                                                       34924, 10732, 35050, 1660,
                                                                       1688, 12412, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37864, 0, 3,
                                                                       35050, 10795, 35176, 1688,
                                                                       1716, 12496, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38032, 0, 3,
                                                                       35176, 10858, 35302, 1716,
                                                                       1744, 12580, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38200, 0, 3,
                                                                       35428, 11110, 35554, 1800,
                                                                       1828, 12832, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38368, 0, 3,
                                                                       35554, 11173, 35680, 1828,
                                                                       1856, 12916, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38536, 0, 3,
                                                                       35680, 11236, 35806, 1856,
                                                                       1884, 13000, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38704, 0, 3,
                                                                       35806, 11299, 35932, 1884,
                                                                       1912, 13084, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38872, 0, 3,
                                                                       35932, 11362, 36058, 1912,
                                                                       1940, 13168, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 39040, 0, 3,
                                                                       36058, 11425, 36184, 1940,
                                                                       1968, 13252, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 39208, 0, 3,
                                                                       36184, 11488, 36310, 1968,
                                                                       1996, 13336, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 39376, 0, 3,
                                                                       36310, 11551, 36436, 1996,
                                                                       2024, 13420, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 39544, 0, 3,
                                                                       36436, 11614, 36562, 2024,
                                                                       2052, 13504, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 39712, 0, 3,
                                                                       36688, 11908, 36856, 2108,
                                                                       2144, 13804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 39928, 0, 3,
                                                                       36856, 11992, 37024, 2144,
                                                                       2180, 13912, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 40144, 0, 3,
                                                                       37024, 12076, 37192, 2180,
                                                                       2216, 14020, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 40360, 0, 3,
                                                                       37192, 12160, 37360, 2216,
                                                                       2252, 14128, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 40576, 0, 3,
                                                                       37360, 12244, 37528, 2252,
                                                                       2288, 14236, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 40792, 0, 3,
                                                                       37528, 12328, 37696, 2288,
                                                                       2324, 14344, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 41008, 0, 3,
                                                                       37696, 12412, 37864, 2324,
                                                                       2360, 14452, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 41224, 0, 3,
                                                                       37864, 12496, 38032, 2360,
                                                                       2396, 14560, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 41440, 0, 3,
                                                                       38200, 12832, 38368, 2468,
                                                                       2504, 14884, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 41656, 0, 3,
                                                                       38368, 12916, 38536, 2504,
                                                                       2540, 14992, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 41872, 0, 3,
                                                                       38536, 13000, 38704, 2540,
                                                                       2576, 15100, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 42088, 0, 3,
                                                                       38704, 13084, 38872, 2576,
                                                                       2612, 15208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 42304, 0, 3,
                                                                       38872, 13168, 39040, 2612,
                                                                       2648, 15316, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 42520, 0, 3,
                                                                       39040, 13252, 39208, 2648,
                                                                       2684, 15424, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 42736, 0, 3,
                                                                       39208, 13336, 39376, 2684,
                                                                       2720, 15532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 42952, 0, 3,
                                                                       39376, 13420, 39544, 2720,
                                                                       2756, 15640, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 43168, 0, 3,
                                                                       39712, 13804, 39928, 2828,
                                                                       2873, 16018, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 43438, 0, 3,
                                                                       39928, 13912, 40144, 2873,
                                                                       2918, 16153, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 43708, 0, 3,
                                                                       40144, 14020, 40360, 2918,
                                                                       2963, 16288, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 43978, 0, 3,
                                                                       40360, 14128, 40576, 2963,
                                                                       3008, 16423, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 44248, 0, 3,
                                                                       40576, 14236, 40792, 3008,
                                                                       3053, 16558, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 44518, 0, 3,
                                                                       40792, 14344, 41008, 3053,
                                                                       3098, 16693, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 44788, 0, 3,
                                                                       41008, 14452, 41224, 3098,
                                                                       3143, 16828, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 45058, 0, 3,
                                                                       41440, 14884, 41656, 3233,
                                                                       3278, 17233, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 45328, 0, 3,
                                                                       41656, 14992, 41872, 3278,
                                                                       3323, 17368, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 45598, 0, 3,
                                                                       41872, 15100, 42088, 3323,
                                                                       3368, 17503, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 45868, 0, 3,
                                                                       42088, 15208, 42304, 3368,
                                                                       3413, 17638, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 46138, 0, 3,
                                                                       42304, 15316, 42520, 3413,
                                                                       3458, 17773, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 46408, 0, 3,
                                                                       42520, 15424, 42736, 3458,
                                                                       3503, 17908, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 46678, 0, 3,
                                                                       42736, 15532, 42952, 3503,
                                                                       3548, 18043, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 46948, 0, 3,
                                                                       43168, 16018, 43438, 3638,
                                                                       3693, 18508, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 47278, 0, 3,
                                                                       43438, 16153, 43708, 3693,
                                                                       3748, 18673, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 47608, 0, 3,
                                                                       43708, 16288, 43978, 3748,
                                                                       3803, 18838, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 47938, 0, 3,
                                                                       43978, 16423, 44248, 3803,
                                                                       3858, 19003, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 48268, 0, 3,
                                                                       44248, 16558, 44518, 3858,
                                                                       3913, 19168, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 48598, 0, 3,
                                                                       44518, 16693, 44788, 3913,
                                                                       3968, 19333, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 48928, 0, 3,
                                                                       45058, 17233, 45328, 4078,
                                                                       4133, 19828, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 49258, 0, 3,
                                                                       45328, 17368, 45598, 4133,
                                                                       4188, 19993, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 49588, 0, 3,
                                                                       45598, 17503, 45868, 4188,
                                                                       4243, 20158, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 49918, 0, 3,
                                                                       45868, 17638, 46138, 4243,
                                                                       4298, 20323, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 50248, 0, 3,
                                                                       46138, 17773, 46408, 4298,
                                                                       4353, 20488, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 50578, 0, 3,
                                                                       46408, 17908, 46678, 4353,
                                                                       4408, 20653, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 50908, 0, 3,
                                                                       46948, 18508, 47278, 4518,
                                                                       4584, 21214, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 51304, 0, 3,
                                                                       47278, 18673, 47608, 4584,
                                                                       4650, 21412, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 51700, 0, 3,
                                                                       47608, 18838, 47938, 4650,
                                                                       4716, 21610, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 52096, 0, 3,
                                                                       47938, 19003, 48268, 4716,
                                                                       4782, 21808, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 52492, 0, 3,
                                                                       48268, 19168, 48598, 4782,
                                                                       4848, 22006, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 52888, 0, 3,
                                                                       48928, 19828, 49258, 4980,
                                                                       5046, 22600, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 53284, 0, 3,
                                                                       49258, 19993, 49588, 5046,
                                                                       5112, 22798, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 53680, 0, 3,
                                                                       49588, 20158, 49918, 5112,
                                                                       5178, 22996, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 54076, 0, 3,
                                                                       49918, 20323, 50248, 5178,
                                                                       5244, 23194, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 54472, 0, 3,
                                                                       50248, 20488, 50578, 5244,
                                                                       5310, 23392, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 54868, 0, 3,
                                                                       50908, 21214, 51304, 5442,
                                                                       5520, 24058, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 55336, 0, 3,
                                                                       51304, 21412, 51700, 5520,
                                                                       5598, 24292, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 55804, 0, 3,
                                                                       51700, 21610, 52096, 5598,
                                                                       5676, 24526, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 56272, 0, 3,
                                                                       52096, 21808, 52492, 5676,
                                                                       5754, 24760, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 56740, 0, 3,
                                                                       52888, 22600, 53284, 5910,
                                                                       5988, 25462, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 57208, 0, 3,
                                                                       53284, 22798, 53680, 5988,
                                                                       6066, 25696, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 57676, 0, 3,
                                                                       53680, 22996, 54076, 6066,
                                                                       6144, 25930, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 58144, 0, 3,
                                                                       54076, 23194, 54472, 6144,
                                                                       6222, 26164, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 58612, 0, 3,
                                                                       54868, 24058, 55336, 6378,
                                                                       6469, 26944, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 59158, 0, 3,
                                                                       55336, 24292, 55804, 6469,
                                                                       6560, 27217, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 59704, 0, 3,
                                                                       55804, 24526, 56272, 6560,
                                                                       6651, 27490, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 60250, 0, 3,
                                                                       56740, 25462, 57208, 6833,
                                                                       6924, 28309, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 60796, 0, 3,
                                                                       57208, 25696, 57676, 6924,
                                                                       7015, 28582, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 61342, 0, 3,
                                                                       57676, 25930, 58144, 7015,
                                                                       7106, 28855, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61888, 3, 7288,
                                                                       7291, 29128, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61898, 3, 7291,
                                                                       7294, 29134, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61908, 3, 7294,
                                                                       7297, 29140, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61918, 3, 7297,
                                                                       7300, 29146, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61928, 3, 7300,
                                                                       7303, 29152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61938, 3, 7303,
                                                                       7306, 29158, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61948, 3, 7306,
                                                                       7309, 29164, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61958, 3, 7309,
                                                                       7312, 29170, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61968, 3, 7312,
                                                                       7315, 29176, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61978, 3, 7315,
                                                                       7318, 29182, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61988, 3, 7318,
                                                                       7321, 29188, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 61998, 3, 7321,
                                                                       7324, 29194, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62008, 3, 7324,
                                                                       7327, 29200, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62018, 3, 7327,
                                                                       7330, 29206, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62028, 3, 7330,
                                                                       7333, 29212, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62038, 3, 7339,
                                                                       7342, 29218, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62048, 3, 7342,
                                                                       7345, 29224, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62058, 3, 7345,
                                                                       7348, 29230, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62068, 3, 7348,
                                                                       7351, 29236, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62078, 3, 7351,
                                                                       7354, 29242, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62088, 3, 7354,
                                                                       7357, 29248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62098, 3, 7357,
                                                                       7360, 29254, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62108, 3, 7360,
                                                                       7363, 29260, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62118, 3, 7363,
                                                                       7366, 29266, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62128, 3, 7366,
                                                                       7369, 29272, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62138, 3, 7369,
                                                                       7372, 29278, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62148, 3, 7372,
                                                                       7375, 29284, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62158, 3, 7375,
                                                                       7378, 29290, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62168, 3, 7378,
                                                                       7381, 29296, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62178, 3, 7381,
                                                                       7384, 29302, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62188, 0, 3,
                                                                       61888, 29128, 61898, 7390,
                                                                       7399, 29308, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62218, 0, 3,
                                                                       61898, 29134, 61908, 7399,
                                                                       7408, 29326, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62248, 0, 3,
                                                                       61908, 29140, 61918, 7408,
                                                                       7417, 29344, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62278, 0, 3,
                                                                       61918, 29146, 61928, 7417,
                                                                       7426, 29362, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62308, 0, 3,
                                                                       61928, 29152, 61938, 7426,
                                                                       7435, 29380, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62338, 0, 3,
                                                                       61938, 29158, 61948, 7435,
                                                                       7444, 29398, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62368, 0, 3,
                                                                       61948, 29164, 61958, 7444,
                                                                       7453, 29416, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62398, 0, 3,
                                                                       61958, 29170, 61968, 7453,
                                                                       7462, 29434, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62428, 0, 3,
                                                                       61968, 29176, 61978, 7462,
                                                                       7471, 29452, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62458, 0, 3,
                                                                       61978, 29182, 61988, 7471,
                                                                       7480, 29470, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62488, 0, 3,
                                                                       61988, 29188, 61998, 7480,
                                                                       7489, 29488, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62518, 0, 3,
                                                                       61998, 29194, 62008, 7489,
                                                                       7498, 29506, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62548, 0, 3,
                                                                       62008, 29200, 62018, 7498,
                                                                       7507, 29524, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62578, 0, 3,
                                                                       62018, 29206, 62028, 7507,
                                                                       7516, 29542, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62608, 0, 3,
                                                                       62038, 29218, 62048, 7534,
                                                                       7543, 29560, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62638, 0, 3,
                                                                       62048, 29224, 62058, 7543,
                                                                       7552, 29578, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62668, 0, 3,
                                                                       62058, 29230, 62068, 7552,
                                                                       7561, 29596, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62698, 0, 3,
                                                                       62068, 29236, 62078, 7561,
                                                                       7570, 29614, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62728, 0, 3,
                                                                       62078, 29242, 62088, 7570,
                                                                       7579, 29632, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62758, 0, 3,
                                                                       62088, 29248, 62098, 7579,
                                                                       7588, 29650, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62788, 0, 3,
                                                                       62098, 29254, 62108, 7588,
                                                                       7597, 29668, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62818, 0, 3,
                                                                       62108, 29260, 62118, 7597,
                                                                       7606, 29686, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62848, 0, 3,
                                                                       62118, 29266, 62128, 7606,
                                                                       7615, 29704, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62878, 0, 3,
                                                                       62128, 29272, 62138, 7615,
                                                                       7624, 29722, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62908, 0, 3,
                                                                       62138, 29278, 62148, 7624,
                                                                       7633, 29740, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62938, 0, 3,
                                                                       62148, 29284, 62158, 7633,
                                                                       7642, 29758, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62968, 0, 3,
                                                                       62158, 29290, 62168, 7642,
                                                                       7651, 29776, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62998, 0, 3,
                                                                       62168, 29296, 62178, 7651,
                                                                       7660, 29794, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63028, 0, 3,
                                                                       62188, 29308, 62218, 7678,
                                                                       7696, 29812, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63088, 0, 3,
                                                                       62218, 29326, 62248, 7696,
                                                                       7714, 29848, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63148, 0, 3,
                                                                       62248, 29344, 62278, 7714,
                                                                       7732, 29884, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63208, 0, 3,
                                                                       62278, 29362, 62308, 7732,
                                                                       7750, 29920, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63268, 0, 3,
                                                                       62308, 29380, 62338, 7750,
                                                                       7768, 29956, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63328, 0, 3,
                                                                       62338, 29398, 62368, 7768,
                                                                       7786, 29992, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63388, 0, 3,
                                                                       62368, 29416, 62398, 7786,
                                                                       7804, 30028, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63448, 0, 3,
                                                                       62398, 29434, 62428, 7804,
                                                                       7822, 30064, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63508, 0, 3,
                                                                       62428, 29452, 62458, 7822,
                                                                       7840, 30100, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63568, 0, 3,
                                                                       62458, 29470, 62488, 7840,
                                                                       7858, 30136, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63628, 0, 3,
                                                                       62488, 29488, 62518, 7858,
                                                                       7876, 30172, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63688, 0, 3,
                                                                       62518, 29506, 62548, 7876,
                                                                       7894, 30208, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63748, 0, 3,
                                                                       62548, 29524, 62578, 7894,
                                                                       7912, 30244, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63808, 0, 3,
                                                                       62608, 29560, 62638, 7948,
                                                                       7966, 30280, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63868, 0, 3,
                                                                       62638, 29578, 62668, 7966,
                                                                       7984, 30316, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63928, 0, 3,
                                                                       62668, 29596, 62698, 7984,
                                                                       8002, 30352, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63988, 0, 3,
                                                                       62698, 29614, 62728, 8002,
                                                                       8020, 30388, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64048, 0, 3,
                                                                       62728, 29632, 62758, 8020,
                                                                       8038, 30424, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64108, 0, 3,
                                                                       62758, 29650, 62788, 8038,
                                                                       8056, 30460, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64168, 0, 3,
                                                                       62788, 29668, 62818, 8056,
                                                                       8074, 30496, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64228, 0, 3,
                                                                       62818, 29686, 62848, 8074,
                                                                       8092, 30532, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64288, 0, 3,
                                                                       62848, 29704, 62878, 8092,
                                                                       8110, 30568, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64348, 0, 3,
                                                                       62878, 29722, 62908, 8110,
                                                                       8128, 30604, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64408, 0, 3,
                                                                       62908, 29740, 62938, 8128,
                                                                       8146, 30640, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64468, 0, 3,
                                                                       62938, 29758, 62968, 8146,
                                                                       8164, 30676, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64528, 0, 3,
                                                                       62968, 29776, 62998, 8164,
                                                                       8182, 30712, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 64588, 0, 3,
                                                                       63028, 29812, 63088, 8218,
                                                                       8248, 30748, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 64688, 0, 3,
                                                                       63088, 29848, 63148, 8248,
                                                                       8278, 30808, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 64788, 0, 3,
                                                                       63148, 29884, 63208, 8278,
                                                                       8308, 30868, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 64888, 0, 3,
                                                                       63208, 29920, 63268, 8308,
                                                                       8338, 30928, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 64988, 0, 3,
                                                                       63268, 29956, 63328, 8338,
                                                                       8368, 30988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65088, 0, 3,
                                                                       63328, 29992, 63388, 8368,
                                                                       8398, 31048, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65188, 0, 3,
                                                                       63388, 30028, 63448, 8398,
                                                                       8428, 31108, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65288, 0, 3,
                                                                       63448, 30064, 63508, 8428,
                                                                       8458, 31168, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65388, 0, 3,
                                                                       63508, 30100, 63568, 8458,
                                                                       8488, 31228, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65488, 0, 3,
                                                                       63568, 30136, 63628, 8488,
                                                                       8518, 31288, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65588, 0, 3,
                                                                       63628, 30172, 63688, 8518,
                                                                       8548, 31348, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65688, 0, 3,
                                                                       63688, 30208, 63748, 8548,
                                                                       8578, 31408, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65788, 0, 3,
                                                                       63808, 30280, 63868, 8638,
                                                                       8668, 31468, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65888, 0, 3,
                                                                       63868, 30316, 63928, 8668,
                                                                       8698, 31528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65988, 0, 3,
                                                                       63928, 30352, 63988, 8698,
                                                                       8728, 31588, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66088, 0, 3,
                                                                       63988, 30388, 64048, 8728,
                                                                       8758, 31648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66188, 0, 3,
                                                                       64048, 30424, 64108, 8758,
                                                                       8788, 31708, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66288, 0, 3,
                                                                       64108, 30460, 64168, 8788,
                                                                       8818, 31768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66388, 0, 3,
                                                                       64168, 30496, 64228, 8818,
                                                                       8848, 31828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66488, 0, 3,
                                                                       64228, 30532, 64288, 8848,
                                                                       8878, 31888, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66588, 0, 3,
                                                                       64288, 30568, 64348, 8878,
                                                                       8908, 31948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66688, 0, 3,
                                                                       64348, 30604, 64408, 8908,
                                                                       8938, 32008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66788, 0, 3,
                                                                       64408, 30640, 64468, 8938,
                                                                       8968, 32068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66888, 0, 3,
                                                                       64468, 30676, 64528, 8968,
                                                                       8998, 32128, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 66988, 0, 3,
                                                                       64588, 30748, 64688, 9058,
                                                                       9103, 32188, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 67138, 0, 3,
                                                                       64688, 30808, 64788, 9103,
                                                                       9148, 32278, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 67288, 0, 3,
                                                                       64788, 30868, 64888, 9148,
                                                                       9193, 32368, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 67438, 0, 3,
                                                                       64888, 30928, 64988, 9193,
                                                                       9238, 32458, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 67588, 0, 3,
                                                                       64988, 30988, 65088, 9238,
                                                                       9283, 32548, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 67738, 0, 3,
                                                                       65088, 31048, 65188, 9283,
                                                                       9328, 32638, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 67888, 0, 3,
                                                                       65188, 31108, 65288, 9328,
                                                                       9373, 32728, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68038, 0, 3,
                                                                       65288, 31168, 65388, 9373,
                                                                       9418, 32818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68188, 0, 3,
                                                                       65388, 31228, 65488, 9418,
                                                                       9463, 32908, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68338, 0, 3,
                                                                       65488, 31288, 65588, 9463,
                                                                       9508, 32998, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68488, 0, 3,
                                                                       65588, 31348, 65688, 9508,
                                                                       9553, 33088, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68638, 0, 3,
                                                                       65788, 31468, 65888, 9643,
                                                                       9688, 33178, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68788, 0, 3,
                                                                       65888, 31528, 65988, 9688,
                                                                       9733, 33268, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68938, 0, 3,
                                                                       65988, 31588, 66088, 9733,
                                                                       9778, 33358, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69088, 0, 3,
                                                                       66088, 31648, 66188, 9778,
                                                                       9823, 33448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69238, 0, 3,
                                                                       66188, 31708, 66288, 9823,
                                                                       9868, 33538, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69388, 0, 3,
                                                                       66288, 31768, 66388, 9868,
                                                                       9913, 33628, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69538, 0, 3,
                                                                       66388, 31828, 66488, 9913,
                                                                       9958, 33718, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69688, 0, 3,
                                                                       66488, 31888, 66588, 9958,
                                                                       10003, 33808, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69838, 0, 3,
                                                                       66588, 31948, 66688,
                                                                       10003, 10048, 33898,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69988, 0, 3,
                                                                       66688, 32008, 66788,
                                                                       10048, 10093, 33988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 70138, 0, 3,
                                                                       66788, 32068, 66888,
                                                                       10093, 10138, 34078,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 70288, 0, 3,
                                                                       66988, 32188, 67138,
                                                                       10228, 10291, 34168,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 70498, 0, 3,
                                                                       67138, 32278, 67288,
                                                                       10291, 10354, 34294,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 70708, 0, 3,
                                                                       67288, 32368, 67438,
                                                                       10354, 10417, 34420,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 70918, 0, 3,
                                                                       67438, 32458, 67588,
                                                                       10417, 10480, 34546,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 71128, 0, 3,
                                                                       67588, 32548, 67738,
                                                                       10480, 10543, 34672,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 71338, 0, 3,
                                                                       67738, 32638, 67888,
                                                                       10543, 10606, 34798,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 71548, 0, 3,
                                                                       67888, 32728, 68038,
                                                                       10606, 10669, 34924,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 71758, 0, 3,
                                                                       68038, 32818, 68188,
                                                                       10669, 10732, 35050,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 71968, 0, 3,
                                                                       68188, 32908, 68338,
                                                                       10732, 10795, 35176,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 72178, 0, 3,
                                                                       68338, 32998, 68488,
                                                                       10795, 10858, 35302,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 72388, 0, 3,
                                                                       68638, 33178, 68788,
                                                                       10984, 11047, 35428,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 72598, 0, 3,
                                                                       68788, 33268, 68938,
                                                                       11047, 11110, 35554,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 72808, 0, 3,
                                                                       68938, 33358, 69088,
                                                                       11110, 11173, 35680,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 73018, 0, 3,
                                                                       69088, 33448, 69238,
                                                                       11173, 11236, 35806,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 73228, 0, 3,
                                                                       69238, 33538, 69388,
                                                                       11236, 11299, 35932,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 73438, 0, 3,
                                                                       69388, 33628, 69538,
                                                                       11299, 11362, 36058,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 73648, 0, 3,
                                                                       69538, 33718, 69688,
                                                                       11362, 11425, 36184,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 73858, 0, 3,
                                                                       69688, 33808, 69838,
                                                                       11425, 11488, 36310,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 74068, 0, 3,
                                                                       69838, 33898, 69988,
                                                                       11488, 11551, 36436,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 74278, 0, 3,
                                                                       69988, 33988, 70138,
                                                                       11551, 11614, 36562,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 74488, 0, 3,
                                                                       70288, 34168, 70498,
                                                                       11740, 11824, 36688,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 74768, 0, 3,
                                                                       70498, 34294, 70708,
                                                                       11824, 11908, 36856,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 75048, 0, 3,
                                                                       70708, 34420, 70918,
                                                                       11908, 11992, 37024,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 75328, 0, 3,
                                                                       70918, 34546, 71128,
                                                                       11992, 12076, 37192,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 75608, 0, 3,
                                                                       71128, 34672, 71338,
                                                                       12076, 12160, 37360,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 75888, 0, 3,
                                                                       71338, 34798, 71548,
                                                                       12160, 12244, 37528,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 76168, 0, 3,
                                                                       71548, 34924, 71758,
                                                                       12244, 12328, 37696,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 76448, 0, 3,
                                                                       71758, 35050, 71968,
                                                                       12328, 12412, 37864,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 76728, 0, 3,
                                                                       71968, 35176, 72178,
                                                                       12412, 12496, 38032,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 77008, 0, 3,
                                                                       72388, 35428, 72598,
                                                                       12664, 12748, 38200,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 77288, 0, 3,
                                                                       72598, 35554, 72808,
                                                                       12748, 12832, 38368,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 77568, 0, 3,
                                                                       72808, 35680, 73018,
                                                                       12832, 12916, 38536,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 77848, 0, 3,
                                                                       73018, 35806, 73228,
                                                                       12916, 13000, 38704,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 78128, 0, 3,
                                                                       73228, 35932, 73438,
                                                                       13000, 13084, 38872,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 78408, 0, 3,
                                                                       73438, 36058, 73648,
                                                                       13084, 13168, 39040,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 78688, 0, 3,
                                                                       73648, 36184, 73858,
                                                                       13168, 13252, 39208,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 78968, 0, 3,
                                                                       73858, 36310, 74068,
                                                                       13252, 13336, 39376,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 79248, 0, 3,
                                                                       74068, 36436, 74278,
                                                                       13336, 13420, 39544,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 79528, 0, 3,
                                                                       74488, 36688, 74768,
                                                                       13588, 13696, 39712,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 79888, 0, 3,
                                                                       74768, 36856, 75048,
                                                                       13696, 13804, 39928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 80248, 0, 3,
                                                                       75048, 37024, 75328,
                                                                       13804, 13912, 40144,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 80608, 0, 3,
                                                                       75328, 37192, 75608,
                                                                       13912, 14020, 40360,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 80968, 0, 3,
                                                                       75608, 37360, 75888,
                                                                       14020, 14128, 40576,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 81328, 0, 3,
                                                                       75888, 37528, 76168,
                                                                       14128, 14236, 40792,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 81688, 0, 3,
                                                                       76168, 37696, 76448,
                                                                       14236, 14344, 41008,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 82048, 0, 3,
                                                                       76448, 37864, 76728,
                                                                       14344, 14452, 41224,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 82408, 0, 3,
                                                                       77008, 38200, 77288,
                                                                       14668, 14776, 41440,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 82768, 0, 3,
                                                                       77288, 38368, 77568,
                                                                       14776, 14884, 41656,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 83128, 0, 3,
                                                                       77568, 38536, 77848,
                                                                       14884, 14992, 41872,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 83488, 0, 3,
                                                                       77848, 38704, 78128,
                                                                       14992, 15100, 42088,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 83848, 0, 3,
                                                                       78128, 38872, 78408,
                                                                       15100, 15208, 42304,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 84208, 0, 3,
                                                                       78408, 39040, 78688,
                                                                       15208, 15316, 42520,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 84568, 0, 3,
                                                                       78688, 39208, 78968,
                                                                       15316, 15424, 42736,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 84928, 0, 3,
                                                                       78968, 39376, 79248,
                                                                       15424, 15532, 42952,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 85288, 0, 3,
                                                                       79528, 39712, 79888,
                                                                       15748, 15883, 43168,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 85738, 0, 3,
                                                                       79888, 39928, 80248,
                                                                       15883, 16018, 43438,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 86188, 0, 3,
                                                                       80248, 40144, 80608,
                                                                       16018, 16153, 43708,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 86638, 0, 3,
                                                                       80608, 40360, 80968,
                                                                       16153, 16288, 43978,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 87088, 0, 3,
                                                                       80968, 40576, 81328,
                                                                       16288, 16423, 44248,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 87538, 0, 3,
                                                                       81328, 40792, 81688,
                                                                       16423, 16558, 44518,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 87988, 0, 3,
                                                                       81688, 41008, 82048,
                                                                       16558, 16693, 44788,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 88438, 0, 3,
                                                                       82408, 41440, 82768,
                                                                       16963, 17098, 45058,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 88888, 0, 3,
                                                                       82768, 41656, 83128,
                                                                       17098, 17233, 45328,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 89338, 0, 3,
                                                                       83128, 41872, 83488,
                                                                       17233, 17368, 45598,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 89788, 0, 3,
                                                                       83488, 42088, 83848,
                                                                       17368, 17503, 45868,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 90238, 0, 3,
                                                                       83848, 42304, 84208,
                                                                       17503, 17638, 46138,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 90688, 0, 3,
                                                                       84208, 42520, 84568,
                                                                       17638, 17773, 46408,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 91138, 0, 3,
                                                                       84568, 42736, 84928,
                                                                       17773, 17908, 46678,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 91588, 0, 3,
                                                                       85288, 43168, 85738,
                                                                       18178, 18343, 46948,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 92138, 0, 3,
                                                                       85738, 43438, 86188,
                                                                       18343, 18508, 47278,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 92688, 0, 3,
                                                                       86188, 43708, 86638,
                                                                       18508, 18673, 47608,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 93238, 0, 3,
                                                                       86638, 43978, 87088,
                                                                       18673, 18838, 47938,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 93788, 0, 3,
                                                                       87088, 44248, 87538,
                                                                       18838, 19003, 48268,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 94338, 0, 3,
                                                                       87538, 44518, 87988,
                                                                       19003, 19168, 48598,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 94888, 0, 3,
                                                                       88438, 45058, 88888,
                                                                       19498, 19663, 48928,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 95438, 0, 3,
                                                                       88888, 45328, 89338,
                                                                       19663, 19828, 49258,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 95988, 0, 3,
                                                                       89338, 45598, 89788,
                                                                       19828, 19993, 49588,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 96538, 0, 3,
                                                                       89788, 45868, 90238,
                                                                       19993, 20158, 49918,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 97088, 0, 3,
                                                                       90238, 46138, 90688,
                                                                       20158, 20323, 50248,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 97638, 0, 3,
                                                                       90688, 46408, 91138,
                                                                       20323, 20488, 50578,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 98188, 0, 3,
                                                                       91588, 46948, 92138,
                                                                       20818, 21016, 50908,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 98848, 0, 3,
                                                                       92138, 47278, 92688,
                                                                       21016, 21214, 51304,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 99508, 0, 3,
                                                                       92688, 47608, 93238,
                                                                       21214, 21412, 51700,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 100168, 0, 3,
                                                                       93238, 47938, 93788,
                                                                       21412, 21610, 52096,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 100828, 0, 3,
                                                                       93788, 48268, 94338,
                                                                       21610, 21808, 52492,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 101488, 0, 3,
                                                                       94888, 48928, 95438,
                                                                       22204, 22402, 52888,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 102148, 0, 3,
                                                                       95438, 49258, 95988,
                                                                       22402, 22600, 53284,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 102808, 0, 3,
                                                                       95988, 49588, 96538,
                                                                       22600, 22798, 53680,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 103468, 0, 3,
                                                                       96538, 49918, 97088,
                                                                       22798, 22996, 54076,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 104128, 0, 3,
                                                                       97088, 50248, 97638,
                                                                       22996, 23194, 54472,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 104788, 0, 3,
                                                                       98188, 50908, 98848,
                                                                       23590, 23824, 54868,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 105568, 0, 3,
                                                                       98848, 51304, 99508,
                                                                       23824, 24058, 55336,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 106348, 0, 3,
                                                                       99508, 51700, 100168,
                                                                       24058, 24292, 55804,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 107128, 0, 3,
                                                                       100168, 52096, 100828,
                                                                       24292, 24526, 56272,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 107908, 0, 3,
                                                                       101488, 52888, 102148,
                                                                       24994, 25228, 56740,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 108688, 0, 3,
                                                                       102148, 53284, 102808,
                                                                       25228, 25462, 57208,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 109468, 0, 3,
                                                                       102808, 53680, 103468,
                                                                       25462, 25696, 57676,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 110248, 0, 3,
                                                                       103468, 54076, 104128,
                                                                       25696, 25930, 58144,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 111028, 0, 3,
                                                                       104788, 54868, 105568,
                                                                       26398, 26671, 58612,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 111938, 0, 3,
                                                                       105568, 55336, 106348,
                                                                       26671, 26944, 59158,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 112848, 0, 3,
                                                                       106348, 55804, 107128,
                                                                       26944, 27217, 59704,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 113758, 0, 3,
                                                                       107908, 56740, 108688,
                                                                       27763, 28036, 60250,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 114668, 0, 3,
                                                                       108688, 57208, 109468,
                                                                       28036, 28309, 60796,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 115578, 0, 3,
                                                                       109468, 57676, 110248,
                                                                       28309, 28582, 61342,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116488, 3, 29128,
                                                                       29134, 61908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116503, 3, 29134,
                                                                       29140, 61918, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116518, 3, 29140,
                                                                       29146, 61928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116533, 3, 29146,
                                                                       29152, 61938, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116548, 3, 29152,
                                                                       29158, 61948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116563, 3, 29158,
                                                                       29164, 61958, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116578, 3, 29164,
                                                                       29170, 61968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116593, 3, 29170,
                                                                       29176, 61978, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116608, 3, 29176,
                                                                       29182, 61988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116623, 3, 29182,
                                                                       29188, 61998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116638, 3, 29188,
                                                                       29194, 62008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116653, 3, 29194,
                                                                       29200, 62018, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116668, 3, 29200,
                                                                       29206, 62028, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116683, 3, 29218,
                                                                       29224, 62058, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116698, 3, 29224,
                                                                       29230, 62068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116713, 3, 29230,
                                                                       29236, 62078, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116728, 3, 29236,
                                                                       29242, 62088, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116743, 3, 29242,
                                                                       29248, 62098, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116758, 3, 29248,
                                                                       29254, 62108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116773, 3, 29254,
                                                                       29260, 62118, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116788, 3, 29260,
                                                                       29266, 62128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116803, 3, 29266,
                                                                       29272, 62138, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116818, 3, 29272,
                                                                       29278, 62148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116833, 3, 29278,
                                                                       29284, 62158, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116848, 3, 29284,
                                                                       29290, 62168, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 116863, 3, 29290,
                                                                       29296, 62178, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 116878, 0, 3,
                                                                       116488, 61908, 116503,
                                                                       29308, 29326, 62248,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 116923, 0, 3,
                                                                       116503, 61918, 116518,
                                                                       29326, 29344, 62278,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 116968, 0, 3,
                                                                       116518, 61928, 116533,
                                                                       29344, 29362, 62308,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117013, 0, 3,
                                                                       116533, 61938, 116548,
                                                                       29362, 29380, 62338,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117058, 0, 3,
                                                                       116548, 61948, 116563,
                                                                       29380, 29398, 62368,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117103, 0, 3,
                                                                       116563, 61958, 116578,
                                                                       29398, 29416, 62398,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117148, 0, 3,
                                                                       116578, 61968, 116593,
                                                                       29416, 29434, 62428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117193, 0, 3,
                                                                       116593, 61978, 116608,
                                                                       29434, 29452, 62458,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117238, 0, 3,
                                                                       116608, 61988, 116623,
                                                                       29452, 29470, 62488,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117283, 0, 3,
                                                                       116623, 61998, 116638,
                                                                       29470, 29488, 62518,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117328, 0, 3,
                                                                       116638, 62008, 116653,
                                                                       29488, 29506, 62548,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117373, 0, 3,
                                                                       116653, 62018, 116668,
                                                                       29506, 29524, 62578,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117418, 0, 3,
                                                                       116683, 62058, 116698,
                                                                       29560, 29578, 62668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117463, 0, 3,
                                                                       116698, 62068, 116713,
                                                                       29578, 29596, 62698,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117508, 0, 3,
                                                                       116713, 62078, 116728,
                                                                       29596, 29614, 62728,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117553, 0, 3,
                                                                       116728, 62088, 116743,
                                                                       29614, 29632, 62758,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117598, 0, 3,
                                                                       116743, 62098, 116758,
                                                                       29632, 29650, 62788,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117643, 0, 3,
                                                                       116758, 62108, 116773,
                                                                       29650, 29668, 62818,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117688, 0, 3,
                                                                       116773, 62118, 116788,
                                                                       29668, 29686, 62848,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117733, 0, 3,
                                                                       116788, 62128, 116803,
                                                                       29686, 29704, 62878,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117778, 0, 3,
                                                                       116803, 62138, 116818,
                                                                       29704, 29722, 62908,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117823, 0, 3,
                                                                       116818, 62148, 116833,
                                                                       29722, 29740, 62938,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117868, 0, 3,
                                                                       116833, 62158, 116848,
                                                                       29740, 29758, 62968,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 117913, 0, 3,
                                                                       116848, 62168, 116863,
                                                                       29758, 29776, 62998,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 117958, 0, 3,
                                                                       116878, 62248, 116923,
                                                                       29812, 29848, 63148,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 118048, 0, 3,
                                                                       116923, 62278, 116968,
                                                                       29848, 29884, 63208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 118138, 0, 3,
                                                                       116968, 62308, 117013,
                                                                       29884, 29920, 63268,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 118228, 0, 3,
                                                                       117013, 62338, 117058,
                                                                       29920, 29956, 63328,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 118318, 0, 3,
                                                                       117058, 62368, 117103,
                                                                       29956, 29992, 63388,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 118408, 0, 3,
                                                                       117103, 62398, 117148,
                                                                       29992, 30028, 63448,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 118498, 0, 3,
                                                                       117148, 62428, 117193,
                                                                       30028, 30064, 63508,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 118588, 0, 3,
                                                                       117193, 62458, 117238,
                                                                       30064, 30100, 63568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 118678, 0, 3,
                                                                       117238, 62488, 117283,
                                                                       30100, 30136, 63628,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 118768, 0, 3,
                                                                       117283, 62518, 117328,
                                                                       30136, 30172, 63688,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 118858, 0, 3,
                                                                       117328, 62548, 117373,
                                                                       30172, 30208, 63748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 118948, 0, 3,
                                                                       117418, 62668, 117463,
                                                                       30280, 30316, 63928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 119038, 0, 3,
                                                                       117463, 62698, 117508,
                                                                       30316, 30352, 63988,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 119128, 0, 3,
                                                                       117508, 62728, 117553,
                                                                       30352, 30388, 64048,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 119218, 0, 3,
                                                                       117553, 62758, 117598,
                                                                       30388, 30424, 64108,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 119308, 0, 3,
                                                                       117598, 62788, 117643,
                                                                       30424, 30460, 64168,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 119398, 0, 3,
                                                                       117643, 62818, 117688,
                                                                       30460, 30496, 64228,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 119488, 0, 3,
                                                                       117688, 62848, 117733,
                                                                       30496, 30532, 64288,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 119578, 0, 3,
                                                                       117733, 62878, 117778,
                                                                       30532, 30568, 64348,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 119668, 0, 3,
                                                                       117778, 62908, 117823,
                                                                       30568, 30604, 64408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 119758, 0, 3,
                                                                       117823, 62938, 117868,
                                                                       30604, 30640, 64468,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 119848, 0, 3,
                                                                       117868, 62968, 117913,
                                                                       30640, 30676, 64528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 119938, 0, 3,
                                                                       117958, 63148, 118048,
                                                                       30748, 30808, 64788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 120088, 0, 3,
                                                                       118048, 63208, 118138,
                                                                       30808, 30868, 64888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 120238, 0, 3,
                                                                       118138, 63268, 118228,
                                                                       30868, 30928, 64988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 120388, 0, 3,
                                                                       118228, 63328, 118318,
                                                                       30928, 30988, 65088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 120538, 0, 3,
                                                                       118318, 63388, 118408,
                                                                       30988, 31048, 65188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 120688, 0, 3,
                                                                       118408, 63448, 118498,
                                                                       31048, 31108, 65288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 120838, 0, 3,
                                                                       118498, 63508, 118588,
                                                                       31108, 31168, 65388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 120988, 0, 3,
                                                                       118588, 63568, 118678,
                                                                       31168, 31228, 65488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 121138, 0, 3,
                                                                       118678, 63628, 118768,
                                                                       31228, 31288, 65588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 121288, 0, 3,
                                                                       118768, 63688, 118858,
                                                                       31288, 31348, 65688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 121438, 0, 3,
                                                                       118948, 63928, 119038,
                                                                       31468, 31528, 65988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 121588, 0, 3,
                                                                       119038, 63988, 119128,
                                                                       31528, 31588, 66088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 121738, 0, 3,
                                                                       119128, 64048, 119218,
                                                                       31588, 31648, 66188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 121888, 0, 3,
                                                                       119218, 64108, 119308,
                                                                       31648, 31708, 66288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 122038, 0, 3,
                                                                       119308, 64168, 119398,
                                                                       31708, 31768, 66388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 122188, 0, 3,
                                                                       119398, 64228, 119488,
                                                                       31768, 31828, 66488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 122338, 0, 3,
                                                                       119488, 64288, 119578,
                                                                       31828, 31888, 66588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 122488, 0, 3,
                                                                       119578, 64348, 119668,
                                                                       31888, 31948, 66688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 122638, 0, 3,
                                                                       119668, 64408, 119758,
                                                                       31948, 32008, 66788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 122788, 0, 3,
                                                                       119758, 64468, 119848,
                                                                       32008, 32068, 66888,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 122938, 0, 3,
                                                                       119938, 64788, 120088,
                                                                       32188, 32278, 67288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 123163, 0, 3,
                                                                       120088, 64888, 120238,
                                                                       32278, 32368, 67438,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 123388, 0, 3,
                                                                       120238, 64988, 120388,
                                                                       32368, 32458, 67588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 123613, 0, 3,
                                                                       120388, 65088, 120538,
                                                                       32458, 32548, 67738,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 123838, 0, 3,
                                                                       120538, 65188, 120688,
                                                                       32548, 32638, 67888,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 124063, 0, 3,
                                                                       120688, 65288, 120838,
                                                                       32638, 32728, 68038,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 124288, 0, 3,
                                                                       120838, 65388, 120988,
                                                                       32728, 32818, 68188,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 124513, 0, 3,
                                                                       120988, 65488, 121138,
                                                                       32818, 32908, 68338,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 124738, 0, 3,
                                                                       121138, 65588, 121288,
                                                                       32908, 32998, 68488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 124963, 0, 3,
                                                                       121438, 65988, 121588,
                                                                       33178, 33268, 68938,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 125188, 0, 3,
                                                                       121588, 66088, 121738,
                                                                       33268, 33358, 69088,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 125413, 0, 3,
                                                                       121738, 66188, 121888,
                                                                       33358, 33448, 69238,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 125638, 0, 3,
                                                                       121888, 66288, 122038,
                                                                       33448, 33538, 69388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 125863, 0, 3,
                                                                       122038, 66388, 122188,
                                                                       33538, 33628, 69538,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 126088, 0, 3,
                                                                       122188, 66488, 122338,
                                                                       33628, 33718, 69688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 126313, 0, 3,
                                                                       122338, 66588, 122488,
                                                                       33718, 33808, 69838,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 126538, 0, 3,
                                                                       122488, 66688, 122638,
                                                                       33808, 33898, 69988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 126763, 0, 3,
                                                                       122638, 66788, 122788,
                                                                       33898, 33988, 70138,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 126988, 0, 3,
                                                                       122938, 67288, 123163,
                                                                       34168, 34294, 70708,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 127303, 0, 3,
                                                                       123163, 67438, 123388,
                                                                       34294, 34420, 70918,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 127618, 0, 3,
                                                                       123388, 67588, 123613,
                                                                       34420, 34546, 71128,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 127933, 0, 3,
                                                                       123613, 67738, 123838,
                                                                       34546, 34672, 71338,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 128248, 0, 3,
                                                                       123838, 67888, 124063,
                                                                       34672, 34798, 71548,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 128563, 0, 3,
                                                                       124063, 68038, 124288,
                                                                       34798, 34924, 71758,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 128878, 0, 3,
                                                                       124288, 68188, 124513,
                                                                       34924, 35050, 71968,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 129193, 0, 3,
                                                                       124513, 68338, 124738,
                                                                       35050, 35176, 72178,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 129508, 0, 3,
                                                                       124963, 68938, 125188,
                                                                       35428, 35554, 72808,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 129823, 0, 3,
                                                                       125188, 69088, 125413,
                                                                       35554, 35680, 73018,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 130138, 0, 3,
                                                                       125413, 69238, 125638,
                                                                       35680, 35806, 73228,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 130453, 0, 3,
                                                                       125638, 69388, 125863,
                                                                       35806, 35932, 73438,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 130768, 0, 3,
                                                                       125863, 69538, 126088,
                                                                       35932, 36058, 73648,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 131083, 0, 3,
                                                                       126088, 69688, 126313,
                                                                       36058, 36184, 73858,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 131398, 0, 3,
                                                                       126313, 69838, 126538,
                                                                       36184, 36310, 74068,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 131713, 0, 3,
                                                                       126538, 69988, 126763,
                                                                       36310, 36436, 74278,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 132028, 0, 3,
                                                                       126988, 70708, 127303,
                                                                       36688, 36856, 75048,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 132448, 0, 3,
                                                                       127303, 70918, 127618,
                                                                       36856, 37024, 75328,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 132868, 0, 3,
                                                                       127618, 71128, 127933,
                                                                       37024, 37192, 75608,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 133288, 0, 3,
                                                                       127933, 71338, 128248,
                                                                       37192, 37360, 75888,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 133708, 0, 3,
                                                                       128248, 71548, 128563,
                                                                       37360, 37528, 76168,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 134128, 0, 3,
                                                                       128563, 71758, 128878,
                                                                       37528, 37696, 76448,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 134548, 0, 3,
                                                                       128878, 71968, 129193,
                                                                       37696, 37864, 76728,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 134968, 0, 3,
                                                                       129508, 72808, 129823,
                                                                       38200, 38368, 77568,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 135388, 0, 3,
                                                                       129823, 73018, 130138,
                                                                       38368, 38536, 77848,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 135808, 0, 3,
                                                                       130138, 73228, 130453,
                                                                       38536, 38704, 78128,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 136228, 0, 3,
                                                                       130453, 73438, 130768,
                                                                       38704, 38872, 78408,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 136648, 0, 3,
                                                                       130768, 73648, 131083,
                                                                       38872, 39040, 78688,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 137068, 0, 3,
                                                                       131083, 73858, 131398,
                                                                       39040, 39208, 78968,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 137488, 0, 3,
                                                                       131398, 74068, 131713,
                                                                       39208, 39376, 79248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 137908, 0, 3,
                                                                       132028, 75048, 132448,
                                                                       39712, 39928, 80248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 138448, 0, 3,
                                                                       132448, 75328, 132868,
                                                                       39928, 40144, 80608,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 138988, 0, 3,
                                                                       132868, 75608, 133288,
                                                                       40144, 40360, 80968,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 139528, 0, 3,
                                                                       133288, 75888, 133708,
                                                                       40360, 40576, 81328,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 140068, 0, 3,
                                                                       133708, 76168, 134128,
                                                                       40576, 40792, 81688,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 140608, 0, 3,
                                                                       134128, 76448, 134548,
                                                                       40792, 41008, 82048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 141148, 0, 3,
                                                                       134968, 77568, 135388,
                                                                       41440, 41656, 83128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 141688, 0, 3,
                                                                       135388, 77848, 135808,
                                                                       41656, 41872, 83488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 142228, 0, 3,
                                                                       135808, 78128, 136228,
                                                                       41872, 42088, 83848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 142768, 0, 3,
                                                                       136228, 78408, 136648,
                                                                       42088, 42304, 84208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 143308, 0, 3,
                                                                       136648, 78688, 137068,
                                                                       42304, 42520, 84568,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 143848, 0, 3,
                                                                       137068, 78968, 137488,
                                                                       42520, 42736, 84928,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 144388, 0, 3,
                                                                       137908, 80248, 138448,
                                                                       43168, 43438, 86188,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 145063, 0, 3,
                                                                       138448, 80608, 138988,
                                                                       43438, 43708, 86638,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 145738, 0, 3,
                                                                       138988, 80968, 139528,
                                                                       43708, 43978, 87088,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 146413, 0, 3,
                                                                       139528, 81328, 140068,
                                                                       43978, 44248, 87538,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 147088, 0, 3,
                                                                       140068, 81688, 140608,
                                                                       44248, 44518, 87988,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 147763, 0, 3,
                                                                       141148, 83128, 141688,
                                                                       45058, 45328, 89338,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 148438, 0, 3,
                                                                       141688, 83488, 142228,
                                                                       45328, 45598, 89788,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 149113, 0, 3,
                                                                       142228, 83848, 142768,
                                                                       45598, 45868, 90238,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 149788, 0, 3,
                                                                       142768, 84208, 143308,
                                                                       45868, 46138, 90688,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 150463, 0, 3,
                                                                       143308, 84568, 143848,
                                                                       46138, 46408, 91138,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 151138, 0, 3,
                                                                       144388, 86188, 145063,
                                                                       46948, 47278, 92688,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 151963, 0, 3,
                                                                       145063, 86638, 145738,
                                                                       47278, 47608, 93238,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 152788, 0, 3,
                                                                       145738, 87088, 146413,
                                                                       47608, 47938, 93788,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 153613, 0, 3,
                                                                       146413, 87538, 147088,
                                                                       47938, 48268, 94338,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 154438, 0, 3,
                                                                       147763, 89338, 148438,
                                                                       48928, 49258, 95988,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 155263, 0, 3,
                                                                       148438, 89788, 149113,
                                                                       49258, 49588, 96538,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 156088, 0, 3,
                                                                       149113, 90238, 149788,
                                                                       49588, 49918, 97088,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 156913, 0, 3,
                                                                       149788, 90688, 150463,
                                                                       49918, 50248, 97638,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 157738, 0, 3,
                                                                       151138, 92688, 151963,
                                                                       50908, 51304, 99508,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 158728, 0, 3,
                                                                       151963, 93238, 152788,
                                                                       51304, 51700, 100168,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 159718, 0, 3,
                                                                       152788, 93788, 153613,
                                                                       51700, 52096, 100828,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 160708, 0, 3,
                                                                       154438, 95988, 155263,
                                                                       52888, 53284, 102808,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 161698, 0, 3,
                                                                       155263, 96538, 156088,
                                                                       53284, 53680, 103468,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 162688, 0, 3,
                                                                       156088, 97088, 156913,
                                                                       53680, 54076, 104128,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 163678, 0, 3,
                                                                       157738, 99508, 158728,
                                                                       54868, 55336, 106348,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 164848, 0, 3,
                                                                       158728, 100168, 159718,
                                                                       55336, 55804, 107128,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 166018, 0, 3,
                                                                       160708, 102808, 161698,
                                                                       56740, 57208, 109468,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 167188, 0, 3,
                                                                       161698, 103468, 162688,
                                                                       57208, 57676, 110248,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 168358, 0, 3,
                                                                       163678, 106348, 164848,
                                                                       58612, 59158, 112848,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 169723, 0, 3,
                                                                       166018, 109468, 167188,
                                                                       60250, 60796, 115578,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171088, 3, 61888,
                                                                       61898, 116488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171109, 3, 61898,
                                                                       61908, 116503, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171130, 3, 61908,
                                                                       61918, 116518, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171151, 3, 61918,
                                                                       61928, 116533, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171172, 3, 61928,
                                                                       61938, 116548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171193, 3, 61938,
                                                                       61948, 116563, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171214, 3, 61948,
                                                                       61958, 116578, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171235, 3, 61958,
                                                                       61968, 116593, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171256, 3, 61968,
                                                                       61978, 116608, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171277, 3, 61978,
                                                                       61988, 116623, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171298, 3, 61988,
                                                                       61998, 116638, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171319, 3, 61998,
                                                                       62008, 116653, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171340, 3, 62008,
                                                                       62018, 116668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171361, 3, 62038,
                                                                       62048, 116683, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171382, 3, 62048,
                                                                       62058, 116698, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171403, 3, 62058,
                                                                       62068, 116713, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171424, 3, 62068,
                                                                       62078, 116728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171445, 3, 62078,
                                                                       62088, 116743, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171466, 3, 62088,
                                                                       62098, 116758, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171487, 3, 62098,
                                                                       62108, 116773, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171508, 3, 62108,
                                                                       62118, 116788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171529, 3, 62118,
                                                                       62128, 116803, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171550, 3, 62128,
                                                                       62138, 116818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171571, 3, 62138,
                                                                       62148, 116833, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171592, 3, 62148,
                                                                       62158, 116848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 171613, 3, 62158,
                                                                       62168, 116863, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 171634, 0, 3,
                                                                       171088, 116488, 171109,
                                                                       62188, 62218, 116878,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 171697, 0, 3,
                                                                       171109, 116503, 171130,
                                                                       62218, 62248, 116923,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 171760, 0, 3,
                                                                       171130, 116518, 171151,
                                                                       62248, 62278, 116968,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 171823, 0, 3,
                                                                       171151, 116533, 171172,
                                                                       62278, 62308, 117013,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 171886, 0, 3,
                                                                       171172, 116548, 171193,
                                                                       62308, 62338, 117058,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 171949, 0, 3,
                                                                       171193, 116563, 171214,
                                                                       62338, 62368, 117103,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172012, 0, 3,
                                                                       171214, 116578, 171235,
                                                                       62368, 62398, 117148,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172075, 0, 3,
                                                                       171235, 116593, 171256,
                                                                       62398, 62428, 117193,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172138, 0, 3,
                                                                       171256, 116608, 171277,
                                                                       62428, 62458, 117238,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172201, 0, 3,
                                                                       171277, 116623, 171298,
                                                                       62458, 62488, 117283,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172264, 0, 3,
                                                                       171298, 116638, 171319,
                                                                       62488, 62518, 117328,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172327, 0, 3,
                                                                       171319, 116653, 171340,
                                                                       62518, 62548, 117373,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172390, 0, 3,
                                                                       171361, 116683, 171382,
                                                                       62608, 62638, 117418,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172453, 0, 3,
                                                                       171382, 116698, 171403,
                                                                       62638, 62668, 117463,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172516, 0, 3,
                                                                       171403, 116713, 171424,
                                                                       62668, 62698, 117508,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172579, 0, 3,
                                                                       171424, 116728, 171445,
                                                                       62698, 62728, 117553,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172642, 0, 3,
                                                                       171445, 116743, 171466,
                                                                       62728, 62758, 117598,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172705, 0, 3,
                                                                       171466, 116758, 171487,
                                                                       62758, 62788, 117643,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172768, 0, 3,
                                                                       171487, 116773, 171508,
                                                                       62788, 62818, 117688,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172831, 0, 3,
                                                                       171508, 116788, 171529,
                                                                       62818, 62848, 117733,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172894, 0, 3,
                                                                       171529, 116803, 171550,
                                                                       62848, 62878, 117778,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 172957, 0, 3,
                                                                       171550, 116818, 171571,
                                                                       62878, 62908, 117823,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 173020, 0, 3,
                                                                       171571, 116833, 171592,
                                                                       62908, 62938, 117868,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 173083, 0, 3,
                                                                       171592, 116848, 171613,
                                                                       62938, 62968, 117913,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 173146, 0, 3,
                                                                       171634, 116878, 171697,
                                                                       63028, 63088, 117958,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 173272, 0, 3,
                                                                       171697, 116923, 171760,
                                                                       63088, 63148, 118048,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 173398, 0, 3,
                                                                       171760, 116968, 171823,
                                                                       63148, 63208, 118138,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 173524, 0, 3,
                                                                       171823, 117013, 171886,
                                                                       63208, 63268, 118228,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 173650, 0, 3,
                                                                       171886, 117058, 171949,
                                                                       63268, 63328, 118318,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 173776, 0, 3,
                                                                       171949, 117103, 172012,
                                                                       63328, 63388, 118408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 173902, 0, 3,
                                                                       172012, 117148, 172075,
                                                                       63388, 63448, 118498,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 174028, 0, 3,
                                                                       172075, 117193, 172138,
                                                                       63448, 63508, 118588,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 174154, 0, 3,
                                                                       172138, 117238, 172201,
                                                                       63508, 63568, 118678,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 174280, 0, 3,
                                                                       172201, 117283, 172264,
                                                                       63568, 63628, 118768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 174406, 0, 3,
                                                                       172264, 117328, 172327,
                                                                       63628, 63688, 118858,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 174532, 0, 3,
                                                                       172390, 117418, 172453,
                                                                       63808, 63868, 118948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 174658, 0, 3,
                                                                       172453, 117463, 172516,
                                                                       63868, 63928, 119038,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 174784, 0, 3,
                                                                       172516, 117508, 172579,
                                                                       63928, 63988, 119128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 174910, 0, 3,
                                                                       172579, 117553, 172642,
                                                                       63988, 64048, 119218,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 175036, 0, 3,
                                                                       172642, 117598, 172705,
                                                                       64048, 64108, 119308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 175162, 0, 3,
                                                                       172705, 117643, 172768,
                                                                       64108, 64168, 119398,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 175288, 0, 3,
                                                                       172768, 117688, 172831,
                                                                       64168, 64228, 119488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 175414, 0, 3,
                                                                       172831, 117733, 172894,
                                                                       64228, 64288, 119578,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 175540, 0, 3,
                                                                       172894, 117778, 172957,
                                                                       64288, 64348, 119668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 175666, 0, 3,
                                                                       172957, 117823, 173020,
                                                                       64348, 64408, 119758,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 175792, 0, 3,
                                                                       173020, 117868, 173083,
                                                                       64408, 64468, 119848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 175918, 0, 3,
                                                                       173146, 117958, 173272,
                                                                       64588, 64688, 119938,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 176128, 0, 3,
                                                                       173272, 118048, 173398,
                                                                       64688, 64788, 120088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 176338, 0, 3,
                                                                       173398, 118138, 173524,
                                                                       64788, 64888, 120238,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 176548, 0, 3,
                                                                       173524, 118228, 173650,
                                                                       64888, 64988, 120388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 176758, 0, 3,
                                                                       173650, 118318, 173776,
                                                                       64988, 65088, 120538,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 176968, 0, 3,
                                                                       173776, 118408, 173902,
                                                                       65088, 65188, 120688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 177178, 0, 3,
                                                                       173902, 118498, 174028,
                                                                       65188, 65288, 120838,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 177388, 0, 3,
                                                                       174028, 118588, 174154,
                                                                       65288, 65388, 120988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 177598, 0, 3,
                                                                       174154, 118678, 174280,
                                                                       65388, 65488, 121138,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 177808, 0, 3,
                                                                       174280, 118768, 174406,
                                                                       65488, 65588, 121288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 178018, 0, 3,
                                                                       174532, 118948, 174658,
                                                                       65788, 65888, 121438,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 178228, 0, 3,
                                                                       174658, 119038, 174784,
                                                                       65888, 65988, 121588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 178438, 0, 3,
                                                                       174784, 119128, 174910,
                                                                       65988, 66088, 121738,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 178648, 0, 3,
                                                                       174910, 119218, 175036,
                                                                       66088, 66188, 121888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 178858, 0, 3,
                                                                       175036, 119308, 175162,
                                                                       66188, 66288, 122038,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 179068, 0, 3,
                                                                       175162, 119398, 175288,
                                                                       66288, 66388, 122188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 179278, 0, 3,
                                                                       175288, 119488, 175414,
                                                                       66388, 66488, 122338,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 179488, 0, 3,
                                                                       175414, 119578, 175540,
                                                                       66488, 66588, 122488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 179698, 0, 3,
                                                                       175540, 119668, 175666,
                                                                       66588, 66688, 122638,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 179908, 0, 3,
                                                                       175666, 119758, 175792,
                                                                       66688, 66788, 122788,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 180118, 0, 3,
                                                                       175918, 119938, 176128,
                                                                       66988, 67138, 122938,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 180433, 0, 3,
                                                                       176128, 120088, 176338,
                                                                       67138, 67288, 123163,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 180748, 0, 3,
                                                                       176338, 120238, 176548,
                                                                       67288, 67438, 123388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 181063, 0, 3,
                                                                       176548, 120388, 176758,
                                                                       67438, 67588, 123613,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 181378, 0, 3,
                                                                       176758, 120538, 176968,
                                                                       67588, 67738, 123838,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 181693, 0, 3,
                                                                       176968, 120688, 177178,
                                                                       67738, 67888, 124063,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 182008, 0, 3,
                                                                       177178, 120838, 177388,
                                                                       67888, 68038, 124288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 182323, 0, 3,
                                                                       177388, 120988, 177598,
                                                                       68038, 68188, 124513,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 182638, 0, 3,
                                                                       177598, 121138, 177808,
                                                                       68188, 68338, 124738,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 182953, 0, 3,
                                                                       178018, 121438, 178228,
                                                                       68638, 68788, 124963,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 183268, 0, 3,
                                                                       178228, 121588, 178438,
                                                                       68788, 68938, 125188,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 183583, 0, 3,
                                                                       178438, 121738, 178648,
                                                                       68938, 69088, 125413,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 183898, 0, 3,
                                                                       178648, 121888, 178858,
                                                                       69088, 69238, 125638,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 184213, 0, 3,
                                                                       178858, 122038, 179068,
                                                                       69238, 69388, 125863,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 184528, 0, 3,
                                                                       179068, 122188, 179278,
                                                                       69388, 69538, 126088,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 184843, 0, 3,
                                                                       179278, 122338, 179488,
                                                                       69538, 69688, 126313,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 185158, 0, 3,
                                                                       179488, 122488, 179698,
                                                                       69688, 69838, 126538,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 185473, 0, 3,
                                                                       179698, 122638, 179908,
                                                                       69838, 69988, 126763,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 185788, 0, 3,
                                                                       180118, 122938, 180433,
                                                                       70288, 70498, 126988,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 186229, 0, 3,
                                                                       180433, 123163, 180748,
                                                                       70498, 70708, 127303,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 186670, 0, 3,
                                                                       180748, 123388, 181063,
                                                                       70708, 70918, 127618,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 187111, 0, 3,
                                                                       181063, 123613, 181378,
                                                                       70918, 71128, 127933,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 187552, 0, 3,
                                                                       181378, 123838, 181693,
                                                                       71128, 71338, 128248,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 187993, 0, 3,
                                                                       181693, 124063, 182008,
                                                                       71338, 71548, 128563,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 188434, 0, 3,
                                                                       182008, 124288, 182323,
                                                                       71548, 71758, 128878,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 188875, 0, 3,
                                                                       182323, 124513, 182638,
                                                                       71758, 71968, 129193,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 189316, 0, 3,
                                                                       182953, 124963, 183268,
                                                                       72388, 72598, 129508,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 189757, 0, 3,
                                                                       183268, 125188, 183583,
                                                                       72598, 72808, 129823,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 190198, 0, 3,
                                                                       183583, 125413, 183898,
                                                                       72808, 73018, 130138,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 190639, 0, 3,
                                                                       183898, 125638, 184213,
                                                                       73018, 73228, 130453,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 191080, 0, 3,
                                                                       184213, 125863, 184528,
                                                                       73228, 73438, 130768,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 191521, 0, 3,
                                                                       184528, 126088, 184843,
                                                                       73438, 73648, 131083,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 191962, 0, 3,
                                                                       184843, 126313, 185158,
                                                                       73648, 73858, 131398,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 192403, 0, 3,
                                                                       185158, 126538, 185473,
                                                                       73858, 74068, 131713,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 192844, 0, 3,
                                                                       185788, 126988, 186229,
                                                                       74488, 74768, 132028,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 193432, 0, 3,
                                                                       186229, 127303, 186670,
                                                                       74768, 75048, 132448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 194020, 0, 3,
                                                                       186670, 127618, 187111,
                                                                       75048, 75328, 132868,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 194608, 0, 3,
                                                                       187111, 127933, 187552,
                                                                       75328, 75608, 133288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 195196, 0, 3,
                                                                       187552, 128248, 187993,
                                                                       75608, 75888, 133708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 195784, 0, 3,
                                                                       187993, 128563, 188434,
                                                                       75888, 76168, 134128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 196372, 0, 3,
                                                                       188434, 128878, 188875,
                                                                       76168, 76448, 134548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 196960, 0, 3,
                                                                       189316, 129508, 189757,
                                                                       77008, 77288, 134968,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 197548, 0, 3,
                                                                       189757, 129823, 190198,
                                                                       77288, 77568, 135388,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 198136, 0, 3,
                                                                       190198, 130138, 190639,
                                                                       77568, 77848, 135808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 198724, 0, 3,
                                                                       190639, 130453, 191080,
                                                                       77848, 78128, 136228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 199312, 0, 3,
                                                                       191080, 130768, 191521,
                                                                       78128, 78408, 136648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 199900, 0, 3,
                                                                       191521, 131083, 191962,
                                                                       78408, 78688, 137068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 200488, 0, 3,
                                                                       191962, 131398, 192403,
                                                                       78688, 78968, 137488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 201076, 0, 3,
                                                                       192844, 132028, 193432,
                                                                       79528, 79888, 137908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 201832, 0, 3,
                                                                       193432, 132448, 194020,
                                                                       79888, 80248, 138448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 202588, 0, 3,
                                                                       194020, 132868, 194608,
                                                                       80248, 80608, 138988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 203344, 0, 3,
                                                                       194608, 133288, 195196,
                                                                       80608, 80968, 139528,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 204100, 0, 3,
                                                                       195196, 133708, 195784,
                                                                       80968, 81328, 140068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 204856, 0, 3,
                                                                       195784, 134128, 196372,
                                                                       81328, 81688, 140608,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 205612, 0, 3,
                                                                       196960, 134968, 197548,
                                                                       82408, 82768, 141148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 206368, 0, 3,
                                                                       197548, 135388, 198136,
                                                                       82768, 83128, 141688,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 207124, 0, 3,
                                                                       198136, 135808, 198724,
                                                                       83128, 83488, 142228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 207880, 0, 3,
                                                                       198724, 136228, 199312,
                                                                       83488, 83848, 142768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 208636, 0, 3,
                                                                       199312, 136648, 199900,
                                                                       83848, 84208, 143308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 209392, 0, 3,
                                                                       199900, 137068, 200488,
                                                                       84208, 84568, 143848,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 210148, 0, 3,
                                                                       201076, 137908, 201832,
                                                                       85288, 85738, 144388,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 211093, 0, 3,
                                                                       201832, 138448, 202588,
                                                                       85738, 86188, 145063,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 212038, 0, 3,
                                                                       202588, 138988, 203344,
                                                                       86188, 86638, 145738,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 212983, 0, 3,
                                                                       203344, 139528, 204100,
                                                                       86638, 87088, 146413,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 213928, 0, 3,
                                                                       204100, 140068, 204856,
                                                                       87088, 87538, 147088,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 214873, 0, 3,
                                                                       205612, 141148, 206368,
                                                                       88438, 88888, 147763,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 215818, 0, 3,
                                                                       206368, 141688, 207124,
                                                                       88888, 89338, 148438,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 216763, 0, 3,
                                                                       207124, 142228, 207880,
                                                                       89338, 89788, 149113,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 217708, 0, 3,
                                                                       207880, 142768, 208636,
                                                                       89788, 90238, 149788,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 218653, 0, 3,
                                                                       208636, 143308, 209392,
                                                                       90238, 90688, 150463,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 219598, 0, 3,
                                                                       210148, 144388, 211093,
                                                                       91588, 92138, 151138,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 220753, 0, 3,
                                                                       211093, 145063, 212038,
                                                                       92138, 92688, 151963,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 221908, 0, 3,
                                                                       212038, 145738, 212983,
                                                                       92688, 93238, 152788,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 223063, 0, 3,
                                                                       212983, 146413, 213928,
                                                                       93238, 93788, 153613,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 224218, 0, 3,
                                                                       214873, 147763, 215818,
                                                                       94888, 95438, 154438,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 225373, 0, 3,
                                                                       215818, 148438, 216763,
                                                                       95438, 95988, 155263,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 226528, 0, 3,
                                                                       216763, 149113, 217708,
                                                                       95988, 96538, 156088,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 227683, 0, 3,
                                                                       217708, 149788, 218653,
                                                                       96538, 97088, 156913,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 228838, 0, 3,
                                                                       219598, 151138, 220753,
                                                                       98188, 98848, 157738,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 230224, 0, 3,
                                                                       220753, 151963, 221908,
                                                                       98848, 99508, 158728,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 231610, 0, 3,
                                                                       221908, 152788, 223063,
                                                                       99508, 100168, 159718,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 232996, 0, 3,
                                                                       224218, 154438, 225373,
                                                                       101488, 102148, 160708,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 234382, 0, 3,
                                                                       225373, 155263, 226528,
                                                                       102148, 102808, 161698,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 235768, 0, 3,
                                                                       226528, 156088, 227683,
                                                                       102808, 103468, 162688,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 237154, 0, 3,
                                                                       228838, 157738, 230224,
                                                                       104788, 105568, 163678,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 238792, 0, 3,
                                                                       230224, 158728, 231610,
                                                                       105568, 106348, 164848,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 240430, 0, 3,
                                                                       232996, 160708, 234382,
                                                                       107908, 108688, 166018,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 242068, 0, 3,
                                                                       234382, 161698, 235768,
                                                                       108688, 109468, 167188,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 243706, 0, 3,
                                                                       237154, 163678, 238792,
                                                                       111028, 111938, 168358,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 245617, 0, 3,
                                                                       240430, 166018, 242068,
                                                                       113758, 114668, 169723,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 247528, 192844, 588, ncols);

                    simdfunc::contract_primitives(buffer, 248424, 196960, 588, ncols);

                    simdfunc::contract_primitives(buffer, 249320, 201076, 756, ncols);

                    simdfunc::contract_primitives(buffer, 250472, 205612, 756, ncols);

                    simdfunc::contract_primitives(buffer, 251624, 210148, 945, ncols);

                    simdfunc::contract_primitives(buffer, 253064, 214873, 945, ncols);

                    simdfunc::contract_primitives(buffer, 254504, 219598, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 256264, 224218, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 258024, 228838, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 260136, 232996, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 262248, 237154, 1638, ncols);

                    simdfunc::contract_primitives(buffer, 264744, 240430, 1638, ncols);

                    simdfunc::contract_primitives(buffer, 267240, 243706, 1911, ncols);

                    simdfunc::contract_primitives(buffer, 270152, 245617, 1911, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 248116, 247528, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 249012, 248424, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 250076, 249320, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 251228, 250472, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 252569, 251624, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 254009, 253064, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 255659, 254504, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 257419, 256264, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 259410, 258024, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 261522, 260136, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 263886, 262248, 78, 1, nmax);

        simdtrf::transform_h_inner(buffer, 266382, 264744, 78, 1, nmax);

        simdtrf::transform_h_inner(buffer, 269151, 267240, 91, 1, nmax);

        simdtrf::transform_h_inner(buffer, 272063, 270152, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 273064, 248116, 250076, 11,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 273988, 249012, 251228, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 274912, 250076, 252569, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 276100, 251228, 254009, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 277288, 252569, 255659, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 278773, 254009, 257419, 11,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 280258, 255659, 259410, 11,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 282073, 257419, 261522, 11,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 283888, 259410, 263886, 11,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 286066, 261522, 266382, 11,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 288244, 263886, 269151, 11,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 290818, 266382, 272063, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 293392, 273064, 274912, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 295240, 273988, 276100, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 297088, 274912, 277288, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 299464, 276100, 278773, 11,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 301840, 277288, 280258, 11,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 304810, 278773, 282073, 11,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 307780, 280258, 283888, 11,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 311410, 282073, 286066, 11,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 315040, 283888, 288244, 11,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 319396, 286066, 290818, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 323752, 293392, 297088, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 326832, 295240, 299464, 11,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 329912, 297088, 301840, 11,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 333872, 299464, 304810, 11,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 337832, 301840, 307780, 11,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 342782, 304810, 311410, 11,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 347732, 307780, 315040, 11,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 353782, 311410, 319396, 11,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 359832, 323752, 329912, 11,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 364452, 326832, 333872, 11,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 369072, 329912, 337832, 11,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 375012, 333872, 342782, 11,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 380952, 337832, 347732, 11,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 388377, 342782, 353782, 11,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 395802, 359832, 369072, 11,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 402270, 364452, 375012, 11,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 408738, 369072, 380952, 11,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 417054, 375012, 388377, 11,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 425370, 395802, 408738, 11,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 433994, 402270, 417054, 11,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 442618, 433994, 28, 11, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 442618, 143, nmax);

        simdtrf::transform_i_inner(buffer, 442618, 425370, 28, 11, nmax);

        simdtrf::transform_i_outer(values + 1859 * nvalues + n * npairs, nvalues, buffer, 442618,
                                   143, nmax);
    }

    for (size_t m = 0; m < 3718; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
