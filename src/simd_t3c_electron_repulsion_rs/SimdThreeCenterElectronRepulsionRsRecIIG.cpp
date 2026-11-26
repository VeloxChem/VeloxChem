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


#include "SimdThreeCenterElectronRepulsionRsRecIIG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_iig_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_iig_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 308582, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3042 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 308582, 147428, 18333, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 16,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 24, 3, 16,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7288, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7291, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7294, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7297, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7300, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7303, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7306, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7309, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7312, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7315, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7318, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7321, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7324, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7327, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7330, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7333, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7336, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7339, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7342, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7345, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7348, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7351, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7354, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7357, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7360, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7363, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7366, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7369, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7372, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7375, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7378, 3, 9, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7387, 3, 10, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7396, 3, 11, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7405, 3, 12, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7414, 3, 13, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7423, 3, 14, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7432, 3, 15, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7441, 3, 16, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7450, 3, 17, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7459, 3, 18, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7468, 3, 19, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7477, 3, 20, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7486, 3, 21, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7495, 3, 22, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7504, 3, 27, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7513, 3, 28, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7522, 3, 29, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7531, 3, 30, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7540, 3, 31, 108,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7549, 3, 32, 111,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7558, 3, 33, 114,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7567, 3, 34, 117,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7576, 3, 35, 120,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7585, 3, 36, 123,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7594, 3, 37, 126,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7603, 3, 38, 129,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7612, 3, 39, 132,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7621, 3, 40, 135,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7630, 3, 48, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7648, 3, 51, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7666, 3, 54, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7684, 3, 57, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7702, 3, 60, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7720, 3, 63, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7738, 3, 66, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7756, 3, 69, 192,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7774, 3, 72, 198,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7792, 3, 75, 204,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7810, 3, 78, 210,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7828, 3, 81, 216,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7846, 3, 84, 222,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7864, 3, 96, 240,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7882, 3, 99, 246,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7900, 3, 102, 252,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7918, 3, 105, 258,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7936, 3, 108, 264,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7954, 3, 111, 270,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7972, 3, 114, 276,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7990, 3, 117, 282,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8008, 3, 120, 288,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8026, 3, 123, 294,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8044, 3, 126, 300,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8062, 3, 129, 306,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8080, 3, 132, 312,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8098, 3, 150, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8128, 3, 156, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8158, 3, 162, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8188, 3, 168, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8218, 3, 174, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8248, 3, 180, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8278, 3, 186, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8308, 3, 192, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8338, 3, 198, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8368, 3, 204, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8398, 3, 210, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8428, 3, 216, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8458, 3, 240, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8488, 3, 246, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8518, 3, 252, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8548, 3, 258, 508,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8578, 3, 264, 518,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8608, 3, 270, 528,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8638, 3, 276, 538,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8668, 3, 282, 548,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8698, 3, 288, 558,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8728, 3, 294, 568,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8758, 3, 300, 578,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8788, 3, 306, 588,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8818, 3, 338, 628,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8863, 3, 348, 643,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8908, 3, 358, 658,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8953, 3, 368, 673,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8998, 3, 378, 688,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9043, 3, 388, 703,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9088, 3, 398, 718,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9133, 3, 408, 733,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9178, 3, 418, 748,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9223, 3, 428, 763,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9268, 3, 438, 778,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9313, 3, 478, 823,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9358, 3, 488, 838,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9403, 3, 498, 853,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9448, 3, 508, 868,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9493, 3, 518, 883,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9538, 3, 528, 898,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9583, 3, 538, 913,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9628, 3, 548, 928,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9673, 3, 558, 943,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9718, 3, 568, 958,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9763, 3, 578, 973,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9808, 3, 628,
                                                                       1030, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9871, 3, 643,
                                                                       1051, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9934, 3, 658,
                                                                       1072, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9997, 3, 673,
                                                                       1093, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10060, 3, 688,
                                                                       1114, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10123, 3, 703,
                                                                       1135, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10186, 3, 718,
                                                                       1156, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10249, 3, 733,
                                                                       1177, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10312, 3, 748,
                                                                       1198, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10375, 3, 763,
                                                                       1219, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10438, 3, 823,
                                                                       1282, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10501, 3, 838,
                                                                       1303, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10564, 3, 853,
                                                                       1324, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10627, 3, 868,
                                                                       1345, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10690, 3, 883,
                                                                       1366, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10753, 3, 898,
                                                                       1387, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10816, 3, 913,
                                                                       1408, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10879, 3, 928,
                                                                       1429, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10942, 3, 943,
                                                                       1450, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11005, 3, 958,
                                                                       1471, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11068, 3, 1030,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11152, 3, 1051,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11236, 3, 1072,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11320, 3, 1093,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11404, 3, 1114,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11488, 3, 1135,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11572, 3, 1156,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11656, 3, 1177,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11740, 3, 1198,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11824, 3, 1282,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11908, 3, 1303,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11992, 3, 1324,
                                                                       1912, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12076, 3, 1345,
                                                                       1940, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12160, 3, 1366,
                                                                       1968, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12244, 3, 1387,
                                                                       1996, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12328, 3, 1408,
                                                                       2024, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12412, 3, 1429,
                                                                       2052, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12496, 3, 1450,
                                                                       2080, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12580, 3, 1548,
                                                                       2180, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12688, 3, 1576,
                                                                       2216, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12796, 3, 1604,
                                                                       2252, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12904, 3, 1632,
                                                                       2288, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13012, 3, 1660,
                                                                       2324, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13120, 3, 1688,
                                                                       2360, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13228, 3, 1716,
                                                                       2396, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13336, 3, 1744,
                                                                       2432, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13444, 3, 1856,
                                                                       2540, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13552, 3, 1884,
                                                                       2576, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13660, 3, 1912,
                                                                       2612, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13768, 3, 1940,
                                                                       2648, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13876, 3, 1968,
                                                                       2684, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13984, 3, 1996,
                                                                       2720, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14092, 3, 2024,
                                                                       2756, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14200, 3, 2052,
                                                                       2792, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14308, 3, 2180,
                                                                       2918, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14443, 3, 2216,
                                                                       2963, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14578, 3, 2252,
                                                                       3008, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14713, 3, 2288,
                                                                       3053, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14848, 3, 2324,
                                                                       3098, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14983, 3, 2360,
                                                                       3143, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15118, 3, 2396,
                                                                       3188, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15253, 3, 2540,
                                                                       3323, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15388, 3, 2576,
                                                                       3368, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15523, 3, 2612,
                                                                       3413, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15658, 3, 2648,
                                                                       3458, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15793, 3, 2684,
                                                                       3503, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15928, 3, 2720,
                                                                       3548, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16063, 3, 2756,
                                                                       3593, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16198, 3, 2918,
                                                                       3748, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16363, 3, 2963,
                                                                       3803, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16528, 3, 3008,
                                                                       3858, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16693, 3, 3053,
                                                                       3913, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16858, 3, 3098,
                                                                       3968, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17023, 3, 3143,
                                                                       4023, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17188, 3, 3323,
                                                                       4188, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17353, 3, 3368,
                                                                       4243, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17518, 3, 3413,
                                                                       4298, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17683, 3, 3458,
                                                                       4353, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17848, 3, 3503,
                                                                       4408, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18013, 3, 3548,
                                                                       4463, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18178, 3, 3748,
                                                                       4650, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18376, 3, 3803,
                                                                       4716, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18574, 3, 3858,
                                                                       4782, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18772, 3, 3913,
                                                                       4848, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18970, 3, 3968,
                                                                       4914, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19168, 3, 4188,
                                                                       5112, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19366, 3, 4243,
                                                                       5178, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19564, 3, 4298,
                                                                       5244, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19762, 3, 4353,
                                                                       5310, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19960, 3, 4408,
                                                                       5376, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 20158, 3, 4650,
                                                                       5598, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 20392, 3, 4716,
                                                                       5676, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 20626, 3, 4782,
                                                                       5754, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 20860, 3, 4848,
                                                                       5832, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 21094, 3, 5112,
                                                                       6066, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 21328, 3, 5178,
                                                                       6144, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 21562, 3, 5244,
                                                                       6222, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 21796, 3, 5310,
                                                                       6300, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 22030, 3, 5598,
                                                                       6560, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 22303, 3, 5676,
                                                                       6651, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 22576, 3, 5754,
                                                                       6742, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 22849, 3, 6066,
                                                                       7015, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 23122, 3, 6144,
                                                                       7106, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 23395, 3, 6222,
                                                                       7197, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23668, 3, 7, 8,
                                                                       7288, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23674, 3, 8, 9,
                                                                       7291, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23680, 3, 9, 10,
                                                                       7294, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23686, 3, 10, 11,
                                                                       7297, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23692, 3, 11, 12,
                                                                       7300, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23698, 3, 12, 13,
                                                                       7303, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23704, 3, 13, 14,
                                                                       7306, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23710, 3, 14, 15,
                                                                       7309, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23716, 3, 15, 16,
                                                                       7312, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23722, 3, 16, 17,
                                                                       7315, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23728, 3, 17, 18,
                                                                       7318, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23734, 3, 18, 19,
                                                                       7321, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23740, 3, 19, 20,
                                                                       7324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23746, 3, 20, 21,
                                                                       7327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23752, 3, 21, 22,
                                                                       7330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23758, 3, 25, 26,
                                                                       7333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23764, 3, 26, 27,
                                                                       7336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23770, 3, 27, 28,
                                                                       7339, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23776, 3, 28, 29,
                                                                       7342, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23782, 3, 29, 30,
                                                                       7345, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23788, 3, 30, 31,
                                                                       7348, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23794, 3, 31, 32,
                                                                       7351, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23800, 3, 32, 33,
                                                                       7354, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23806, 3, 33, 34,
                                                                       7357, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23812, 3, 34, 35,
                                                                       7360, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23818, 3, 35, 36,
                                                                       7363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23824, 3, 36, 37,
                                                                       7366, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23830, 3, 37, 38,
                                                                       7369, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23836, 3, 38, 39,
                                                                       7372, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 23842, 3, 39, 40,
                                                                       7375, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23848, 0, 3,
                                                                       23668, 7288, 23674, 7378,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23866, 0, 3,
                                                                       23674, 7291, 23680, 7387,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23884, 0, 3,
                                                                       23680, 7294, 23686, 7396,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23902, 0, 3,
                                                                       23686, 7297, 23692, 7405,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23920, 0, 3,
                                                                       23692, 7300, 23698, 7414,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23938, 0, 3,
                                                                       23698, 7303, 23704, 7423,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23956, 0, 3,
                                                                       23704, 7306, 23710, 7432,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23974, 0, 3,
                                                                       23710, 7309, 23716, 7441,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23992, 0, 3,
                                                                       23716, 7312, 23722, 7450,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24010, 0, 3,
                                                                       23722, 7315, 23728, 7459,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24028, 0, 3,
                                                                       23728, 7318, 23734, 7468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24046, 0, 3,
                                                                       23734, 7321, 23740, 7477,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24064, 0, 3,
                                                                       23740, 7324, 23746, 7486,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24082, 0, 3,
                                                                       23746, 7327, 23752, 7495,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24100, 0, 3,
                                                                       23758, 7333, 23764, 7504,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24118, 0, 3,
                                                                       23764, 7336, 23770, 7513,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24136, 0, 3,
                                                                       23770, 7339, 23776, 7522,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24154, 0, 3,
                                                                       23776, 7342, 23782, 7531,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24172, 0, 3,
                                                                       23782, 7345, 23788, 7540,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24190, 0, 3,
                                                                       23788, 7348, 23794, 7549,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24208, 0, 3,
                                                                       23794, 7351, 23800, 7558,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24226, 0, 3,
                                                                       23800, 7354, 23806, 7567,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24244, 0, 3,
                                                                       23806, 7357, 23812, 7576,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24262, 0, 3,
                                                                       23812, 7360, 23818, 7585,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24280, 0, 3,
                                                                       23818, 7363, 23824, 7594,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24298, 0, 3,
                                                                       23824, 7366, 23830, 7603,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24316, 0, 3,
                                                                       23830, 7369, 23836, 7612,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24334, 0, 3,
                                                                       23836, 7372, 23842, 7621,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24352, 0, 3,
                                                                       23848, 7378, 23866, 138,
                                                                       144, 7630, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24388, 0, 3,
                                                                       23866, 7387, 23884, 144,
                                                                       150, 7648, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24424, 0, 3,
                                                                       23884, 7396, 23902, 150,
                                                                       156, 7666, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24460, 0, 3,
                                                                       23902, 7405, 23920, 156,
                                                                       162, 7684, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24496, 0, 3,
                                                                       23920, 7414, 23938, 162,
                                                                       168, 7702, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24532, 0, 3,
                                                                       23938, 7423, 23956, 168,
                                                                       174, 7720, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24568, 0, 3,
                                                                       23956, 7432, 23974, 174,
                                                                       180, 7738, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24604, 0, 3,
                                                                       23974, 7441, 23992, 180,
                                                                       186, 7756, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24640, 0, 3,
                                                                       23992, 7450, 24010, 186,
                                                                       192, 7774, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24676, 0, 3,
                                                                       24010, 7459, 24028, 192,
                                                                       198, 7792, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24712, 0, 3,
                                                                       24028, 7468, 24046, 198,
                                                                       204, 7810, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24748, 0, 3,
                                                                       24046, 7477, 24064, 204,
                                                                       210, 7828, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24784, 0, 3,
                                                                       24064, 7486, 24082, 210,
                                                                       216, 7846, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24820, 0, 3,
                                                                       24100, 7504, 24118, 228,
                                                                       234, 7864, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24856, 0, 3,
                                                                       24118, 7513, 24136, 234,
                                                                       240, 7882, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24892, 0, 3,
                                                                       24136, 7522, 24154, 240,
                                                                       246, 7900, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24928, 0, 3,
                                                                       24154, 7531, 24172, 246,
                                                                       252, 7918, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24964, 0, 3,
                                                                       24172, 7540, 24190, 252,
                                                                       258, 7936, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25000, 0, 3,
                                                                       24190, 7549, 24208, 258,
                                                                       264, 7954, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25036, 0, 3,
                                                                       24208, 7558, 24226, 264,
                                                                       270, 7972, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25072, 0, 3,
                                                                       24226, 7567, 24244, 270,
                                                                       276, 7990, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25108, 0, 3,
                                                                       24244, 7576, 24262, 276,
                                                                       282, 8008, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25144, 0, 3,
                                                                       24262, 7585, 24280, 282,
                                                                       288, 8026, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25180, 0, 3,
                                                                       24280, 7594, 24298, 288,
                                                                       294, 8044, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25216, 0, 3,
                                                                       24298, 7603, 24316, 294,
                                                                       300, 8062, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25252, 0, 3,
                                                                       24316, 7612, 24334, 300,
                                                                       306, 8080, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25288, 0, 3,
                                                                       24352, 7630, 24388, 318,
                                                                       328, 8098, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25348, 0, 3,
                                                                       24388, 7648, 24424, 328,
                                                                       338, 8128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25408, 0, 3,
                                                                       24424, 7666, 24460, 338,
                                                                       348, 8158, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25468, 0, 3,
                                                                       24460, 7684, 24496, 348,
                                                                       358, 8188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25528, 0, 3,
                                                                       24496, 7702, 24532, 358,
                                                                       368, 8218, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25588, 0, 3,
                                                                       24532, 7720, 24568, 368,
                                                                       378, 8248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25648, 0, 3,
                                                                       24568, 7738, 24604, 378,
                                                                       388, 8278, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25708, 0, 3,
                                                                       24604, 7756, 24640, 388,
                                                                       398, 8308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25768, 0, 3,
                                                                       24640, 7774, 24676, 398,
                                                                       408, 8338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25828, 0, 3,
                                                                       24676, 7792, 24712, 408,
                                                                       418, 8368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25888, 0, 3,
                                                                       24712, 7810, 24748, 418,
                                                                       428, 8398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25948, 0, 3,
                                                                       24748, 7828, 24784, 428,
                                                                       438, 8428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26008, 0, 3,
                                                                       24820, 7864, 24856, 458,
                                                                       468, 8458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26068, 0, 3,
                                                                       24856, 7882, 24892, 468,
                                                                       478, 8488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26128, 0, 3,
                                                                       24892, 7900, 24928, 478,
                                                                       488, 8518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26188, 0, 3,
                                                                       24928, 7918, 24964, 488,
                                                                       498, 8548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26248, 0, 3,
                                                                       24964, 7936, 25000, 498,
                                                                       508, 8578, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26308, 0, 3,
                                                                       25000, 7954, 25036, 508,
                                                                       518, 8608, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26368, 0, 3,
                                                                       25036, 7972, 25072, 518,
                                                                       528, 8638, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26428, 0, 3,
                                                                       25072, 7990, 25108, 528,
                                                                       538, 8668, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26488, 0, 3,
                                                                       25108, 8008, 25144, 538,
                                                                       548, 8698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26548, 0, 3,
                                                                       25144, 8026, 25180, 548,
                                                                       558, 8728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26608, 0, 3,
                                                                       25180, 8044, 25216, 558,
                                                                       568, 8758, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26668, 0, 3,
                                                                       25216, 8062, 25252, 568,
                                                                       578, 8788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26728, 0, 3,
                                                                       25288, 8098, 25348, 598,
                                                                       613, 8818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26818, 0, 3,
                                                                       25348, 8128, 25408, 613,
                                                                       628, 8863, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26908, 0, 3,
                                                                       25408, 8158, 25468, 628,
                                                                       643, 8908, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26998, 0, 3,
                                                                       25468, 8188, 25528, 643,
                                                                       658, 8953, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27088, 0, 3,
                                                                       25528, 8218, 25588, 658,
                                                                       673, 8998, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27178, 0, 3,
                                                                       25588, 8248, 25648, 673,
                                                                       688, 9043, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27268, 0, 3,
                                                                       25648, 8278, 25708, 688,
                                                                       703, 9088, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27358, 0, 3,
                                                                       25708, 8308, 25768, 703,
                                                                       718, 9133, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27448, 0, 3,
                                                                       25768, 8338, 25828, 718,
                                                                       733, 9178, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27538, 0, 3,
                                                                       25828, 8368, 25888, 733,
                                                                       748, 9223, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27628, 0, 3,
                                                                       25888, 8398, 25948, 748,
                                                                       763, 9268, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27718, 0, 3,
                                                                       26008, 8458, 26068, 793,
                                                                       808, 9313, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27808, 0, 3,
                                                                       26068, 8488, 26128, 808,
                                                                       823, 9358, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27898, 0, 3,
                                                                       26128, 8518, 26188, 823,
                                                                       838, 9403, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27988, 0, 3,
                                                                       26188, 8548, 26248, 838,
                                                                       853, 9448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28078, 0, 3,
                                                                       26248, 8578, 26308, 853,
                                                                       868, 9493, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28168, 0, 3,
                                                                       26308, 8608, 26368, 868,
                                                                       883, 9538, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28258, 0, 3,
                                                                       26368, 8638, 26428, 883,
                                                                       898, 9583, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28348, 0, 3,
                                                                       26428, 8668, 26488, 898,
                                                                       913, 9628, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28438, 0, 3,
                                                                       26488, 8698, 26548, 913,
                                                                       928, 9673, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28528, 0, 3,
                                                                       26548, 8728, 26608, 928,
                                                                       943, 9718, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28618, 0, 3,
                                                                       26608, 8758, 26668, 943,
                                                                       958, 9763, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28708, 0, 3,
                                                                       26728, 8818, 26818, 988,
                                                                       1009, 9808, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28834, 0, 3,
                                                                       26818, 8863, 26908, 1009,
                                                                       1030, 9871, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28960, 0, 3,
                                                                       26908, 8908, 26998, 1030,
                                                                       1051, 9934, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29086, 0, 3,
                                                                       26998, 8953, 27088, 1051,
                                                                       1072, 9997, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29212, 0, 3,
                                                                       27088, 8998, 27178, 1072,
                                                                       1093, 10060, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29338, 0, 3,
                                                                       27178, 9043, 27268, 1093,
                                                                       1114, 10123, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29464, 0, 3,
                                                                       27268, 9088, 27358, 1114,
                                                                       1135, 10186, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29590, 0, 3,
                                                                       27358, 9133, 27448, 1135,
                                                                       1156, 10249, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29716, 0, 3,
                                                                       27448, 9178, 27538, 1156,
                                                                       1177, 10312, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29842, 0, 3,
                                                                       27538, 9223, 27628, 1177,
                                                                       1198, 10375, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29968, 0, 3,
                                                                       27718, 9313, 27808, 1240,
                                                                       1261, 10438, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30094, 0, 3,
                                                                       27808, 9358, 27898, 1261,
                                                                       1282, 10501, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30220, 0, 3,
                                                                       27898, 9403, 27988, 1282,
                                                                       1303, 10564, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30346, 0, 3,
                                                                       27988, 9448, 28078, 1303,
                                                                       1324, 10627, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30472, 0, 3,
                                                                       28078, 9493, 28168, 1324,
                                                                       1345, 10690, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30598, 0, 3,
                                                                       28168, 9538, 28258, 1345,
                                                                       1366, 10753, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30724, 0, 3,
                                                                       28258, 9583, 28348, 1366,
                                                                       1387, 10816, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30850, 0, 3,
                                                                       28348, 9628, 28438, 1387,
                                                                       1408, 10879, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30976, 0, 3,
                                                                       28438, 9673, 28528, 1408,
                                                                       1429, 10942, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 31102, 0, 3,
                                                                       28528, 9718, 28618, 1429,
                                                                       1450, 11005, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31228, 0, 3,
                                                                       28708, 9808, 28834, 1492,
                                                                       1520, 11068, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31396, 0, 3,
                                                                       28834, 9871, 28960, 1520,
                                                                       1548, 11152, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31564, 0, 3,
                                                                       28960, 9934, 29086, 1548,
                                                                       1576, 11236, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31732, 0, 3,
                                                                       29086, 9997, 29212, 1576,
                                                                       1604, 11320, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31900, 0, 3,
                                                                       29212, 10060, 29338, 1604,
                                                                       1632, 11404, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32068, 0, 3,
                                                                       29338, 10123, 29464, 1632,
                                                                       1660, 11488, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32236, 0, 3,
                                                                       29464, 10186, 29590, 1660,
                                                                       1688, 11572, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32404, 0, 3,
                                                                       29590, 10249, 29716, 1688,
                                                                       1716, 11656, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32572, 0, 3,
                                                                       29716, 10312, 29842, 1716,
                                                                       1744, 11740, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32740, 0, 3,
                                                                       29968, 10438, 30094, 1800,
                                                                       1828, 11824, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32908, 0, 3,
                                                                       30094, 10501, 30220, 1828,
                                                                       1856, 11908, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33076, 0, 3,
                                                                       30220, 10564, 30346, 1856,
                                                                       1884, 11992, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33244, 0, 3,
                                                                       30346, 10627, 30472, 1884,
                                                                       1912, 12076, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33412, 0, 3,
                                                                       30472, 10690, 30598, 1912,
                                                                       1940, 12160, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33580, 0, 3,
                                                                       30598, 10753, 30724, 1940,
                                                                       1968, 12244, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33748, 0, 3,
                                                                       30724, 10816, 30850, 1968,
                                                                       1996, 12328, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33916, 0, 3,
                                                                       30850, 10879, 30976, 1996,
                                                                       2024, 12412, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 34084, 0, 3,
                                                                       30976, 10942, 31102, 2024,
                                                                       2052, 12496, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34252, 0, 3,
                                                                       31228, 11068, 31396, 2108,
                                                                       2144, 12580, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34468, 0, 3,
                                                                       31396, 11152, 31564, 2144,
                                                                       2180, 12688, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34684, 0, 3,
                                                                       31564, 11236, 31732, 2180,
                                                                       2216, 12796, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34900, 0, 3,
                                                                       31732, 11320, 31900, 2216,
                                                                       2252, 12904, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 35116, 0, 3,
                                                                       31900, 11404, 32068, 2252,
                                                                       2288, 13012, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 35332, 0, 3,
                                                                       32068, 11488, 32236, 2288,
                                                                       2324, 13120, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 35548, 0, 3,
                                                                       32236, 11572, 32404, 2324,
                                                                       2360, 13228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 35764, 0, 3,
                                                                       32404, 11656, 32572, 2360,
                                                                       2396, 13336, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 35980, 0, 3,
                                                                       32740, 11824, 32908, 2468,
                                                                       2504, 13444, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 36196, 0, 3,
                                                                       32908, 11908, 33076, 2504,
                                                                       2540, 13552, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 36412, 0, 3,
                                                                       33076, 11992, 33244, 2540,
                                                                       2576, 13660, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 36628, 0, 3,
                                                                       33244, 12076, 33412, 2576,
                                                                       2612, 13768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 36844, 0, 3,
                                                                       33412, 12160, 33580, 2612,
                                                                       2648, 13876, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 37060, 0, 3,
                                                                       33580, 12244, 33748, 2648,
                                                                       2684, 13984, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 37276, 0, 3,
                                                                       33748, 12328, 33916, 2684,
                                                                       2720, 14092, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 37492, 0, 3,
                                                                       33916, 12412, 34084, 2720,
                                                                       2756, 14200, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 37708, 0, 3,
                                                                       34252, 12580, 34468, 2828,
                                                                       2873, 14308, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 37978, 0, 3,
                                                                       34468, 12688, 34684, 2873,
                                                                       2918, 14443, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 38248, 0, 3,
                                                                       34684, 12796, 34900, 2918,
                                                                       2963, 14578, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 38518, 0, 3,
                                                                       34900, 12904, 35116, 2963,
                                                                       3008, 14713, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 38788, 0, 3,
                                                                       35116, 13012, 35332, 3008,
                                                                       3053, 14848, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 39058, 0, 3,
                                                                       35332, 13120, 35548, 3053,
                                                                       3098, 14983, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 39328, 0, 3,
                                                                       35548, 13228, 35764, 3098,
                                                                       3143, 15118, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 39598, 0, 3,
                                                                       35980, 13444, 36196, 3233,
                                                                       3278, 15253, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 39868, 0, 3,
                                                                       36196, 13552, 36412, 3278,
                                                                       3323, 15388, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 40138, 0, 3,
                                                                       36412, 13660, 36628, 3323,
                                                                       3368, 15523, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 40408, 0, 3,
                                                                       36628, 13768, 36844, 3368,
                                                                       3413, 15658, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 40678, 0, 3,
                                                                       36844, 13876, 37060, 3413,
                                                                       3458, 15793, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 40948, 0, 3,
                                                                       37060, 13984, 37276, 3458,
                                                                       3503, 15928, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 41218, 0, 3,
                                                                       37276, 14092, 37492, 3503,
                                                                       3548, 16063, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 41488, 0, 3,
                                                                       37708, 14308, 37978, 3638,
                                                                       3693, 16198, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 41818, 0, 3,
                                                                       37978, 14443, 38248, 3693,
                                                                       3748, 16363, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 42148, 0, 3,
                                                                       38248, 14578, 38518, 3748,
                                                                       3803, 16528, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 42478, 0, 3,
                                                                       38518, 14713, 38788, 3803,
                                                                       3858, 16693, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 42808, 0, 3,
                                                                       38788, 14848, 39058, 3858,
                                                                       3913, 16858, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 43138, 0, 3,
                                                                       39058, 14983, 39328, 3913,
                                                                       3968, 17023, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 43468, 0, 3,
                                                                       39598, 15253, 39868, 4078,
                                                                       4133, 17188, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 43798, 0, 3,
                                                                       39868, 15388, 40138, 4133,
                                                                       4188, 17353, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 44128, 0, 3,
                                                                       40138, 15523, 40408, 4188,
                                                                       4243, 17518, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 44458, 0, 3,
                                                                       40408, 15658, 40678, 4243,
                                                                       4298, 17683, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 44788, 0, 3,
                                                                       40678, 15793, 40948, 4298,
                                                                       4353, 17848, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 45118, 0, 3,
                                                                       40948, 15928, 41218, 4353,
                                                                       4408, 18013, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 45448, 0, 3,
                                                                       41488, 16198, 41818, 4518,
                                                                       4584, 18178, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 45844, 0, 3,
                                                                       41818, 16363, 42148, 4584,
                                                                       4650, 18376, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 46240, 0, 3,
                                                                       42148, 16528, 42478, 4650,
                                                                       4716, 18574, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 46636, 0, 3,
                                                                       42478, 16693, 42808, 4716,
                                                                       4782, 18772, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 47032, 0, 3,
                                                                       42808, 16858, 43138, 4782,
                                                                       4848, 18970, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 47428, 0, 3,
                                                                       43468, 17188, 43798, 4980,
                                                                       5046, 19168, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 47824, 0, 3,
                                                                       43798, 17353, 44128, 5046,
                                                                       5112, 19366, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 48220, 0, 3,
                                                                       44128, 17518, 44458, 5112,
                                                                       5178, 19564, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 48616, 0, 3,
                                                                       44458, 17683, 44788, 5178,
                                                                       5244, 19762, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 49012, 0, 3,
                                                                       44788, 17848, 45118, 5244,
                                                                       5310, 19960, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 49408, 0, 3,
                                                                       45448, 18178, 45844, 5442,
                                                                       5520, 20158, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 49876, 0, 3,
                                                                       45844, 18376, 46240, 5520,
                                                                       5598, 20392, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 50344, 0, 3,
                                                                       46240, 18574, 46636, 5598,
                                                                       5676, 20626, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 50812, 0, 3,
                                                                       46636, 18772, 47032, 5676,
                                                                       5754, 20860, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 51280, 0, 3,
                                                                       47428, 19168, 47824, 5910,
                                                                       5988, 21094, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 51748, 0, 3,
                                                                       47824, 19366, 48220, 5988,
                                                                       6066, 21328, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 52216, 0, 3,
                                                                       48220, 19564, 48616, 6066,
                                                                       6144, 21562, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 52684, 0, 3,
                                                                       48616, 19762, 49012, 6144,
                                                                       6222, 21796, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 53152, 0, 3,
                                                                       49408, 20158, 49876, 6378,
                                                                       6469, 22030, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 53698, 0, 3,
                                                                       49876, 20392, 50344, 6469,
                                                                       6560, 22303, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 54244, 0, 3,
                                                                       50344, 20626, 50812, 6560,
                                                                       6651, 22576, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 54790, 0, 3,
                                                                       51280, 21094, 51748, 6833,
                                                                       6924, 22849, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 55336, 0, 3,
                                                                       51748, 21328, 52216, 6924,
                                                                       7015, 23122, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 55882, 0, 3,
                                                                       52216, 21562, 52684, 7015,
                                                                       7106, 23395, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56428, 3, 7288,
                                                                       7291, 23680, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56438, 3, 7291,
                                                                       7294, 23686, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56448, 3, 7294,
                                                                       7297, 23692, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56458, 3, 7297,
                                                                       7300, 23698, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56468, 3, 7300,
                                                                       7303, 23704, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56478, 3, 7303,
                                                                       7306, 23710, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56488, 3, 7306,
                                                                       7309, 23716, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56498, 3, 7309,
                                                                       7312, 23722, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56508, 3, 7312,
                                                                       7315, 23728, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56518, 3, 7315,
                                                                       7318, 23734, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56528, 3, 7318,
                                                                       7321, 23740, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56538, 3, 7321,
                                                                       7324, 23746, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56548, 3, 7324,
                                                                       7327, 23752, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56558, 3, 7333,
                                                                       7336, 23770, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56568, 3, 7336,
                                                                       7339, 23776, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56578, 3, 7339,
                                                                       7342, 23782, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56588, 3, 7342,
                                                                       7345, 23788, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56598, 3, 7345,
                                                                       7348, 23794, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56608, 3, 7348,
                                                                       7351, 23800, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56618, 3, 7351,
                                                                       7354, 23806, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56628, 3, 7354,
                                                                       7357, 23812, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56638, 3, 7357,
                                                                       7360, 23818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56648, 3, 7360,
                                                                       7363, 23824, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56658, 3, 7363,
                                                                       7366, 23830, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56668, 3, 7366,
                                                                       7369, 23836, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 56678, 3, 7369,
                                                                       7372, 23842, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 56688, 0, 3,
                                                                       56428, 23680, 56438, 7378,
                                                                       7387, 23884, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 56718, 0, 3,
                                                                       56438, 23686, 56448, 7387,
                                                                       7396, 23902, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 56748, 0, 3,
                                                                       56448, 23692, 56458, 7396,
                                                                       7405, 23920, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 56778, 0, 3,
                                                                       56458, 23698, 56468, 7405,
                                                                       7414, 23938, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 56808, 0, 3,
                                                                       56468, 23704, 56478, 7414,
                                                                       7423, 23956, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 56838, 0, 3,
                                                                       56478, 23710, 56488, 7423,
                                                                       7432, 23974, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 56868, 0, 3,
                                                                       56488, 23716, 56498, 7432,
                                                                       7441, 23992, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 56898, 0, 3,
                                                                       56498, 23722, 56508, 7441,
                                                                       7450, 24010, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 56928, 0, 3,
                                                                       56508, 23728, 56518, 7450,
                                                                       7459, 24028, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 56958, 0, 3,
                                                                       56518, 23734, 56528, 7459,
                                                                       7468, 24046, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 56988, 0, 3,
                                                                       56528, 23740, 56538, 7468,
                                                                       7477, 24064, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57018, 0, 3,
                                                                       56538, 23746, 56548, 7477,
                                                                       7486, 24082, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57048, 0, 3,
                                                                       56558, 23770, 56568, 7504,
                                                                       7513, 24136, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57078, 0, 3,
                                                                       56568, 23776, 56578, 7513,
                                                                       7522, 24154, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57108, 0, 3,
                                                                       56578, 23782, 56588, 7522,
                                                                       7531, 24172, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57138, 0, 3,
                                                                       56588, 23788, 56598, 7531,
                                                                       7540, 24190, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57168, 0, 3,
                                                                       56598, 23794, 56608, 7540,
                                                                       7549, 24208, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57198, 0, 3,
                                                                       56608, 23800, 56618, 7549,
                                                                       7558, 24226, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57228, 0, 3,
                                                                       56618, 23806, 56628, 7558,
                                                                       7567, 24244, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57258, 0, 3,
                                                                       56628, 23812, 56638, 7567,
                                                                       7576, 24262, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57288, 0, 3,
                                                                       56638, 23818, 56648, 7576,
                                                                       7585, 24280, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57318, 0, 3,
                                                                       56648, 23824, 56658, 7585,
                                                                       7594, 24298, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57348, 0, 3,
                                                                       56658, 23830, 56668, 7594,
                                                                       7603, 24316, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 57378, 0, 3,
                                                                       56668, 23836, 56678, 7603,
                                                                       7612, 24334, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 57408, 0, 3,
                                                                       56688, 23884, 56718, 7630,
                                                                       7648, 24424, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 57468, 0, 3,
                                                                       56718, 23902, 56748, 7648,
                                                                       7666, 24460, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 57528, 0, 3,
                                                                       56748, 23920, 56778, 7666,
                                                                       7684, 24496, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 57588, 0, 3,
                                                                       56778, 23938, 56808, 7684,
                                                                       7702, 24532, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 57648, 0, 3,
                                                                       56808, 23956, 56838, 7702,
                                                                       7720, 24568, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 57708, 0, 3,
                                                                       56838, 23974, 56868, 7720,
                                                                       7738, 24604, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 57768, 0, 3,
                                                                       56868, 23992, 56898, 7738,
                                                                       7756, 24640, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 57828, 0, 3,
                                                                       56898, 24010, 56928, 7756,
                                                                       7774, 24676, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 57888, 0, 3,
                                                                       56928, 24028, 56958, 7774,
                                                                       7792, 24712, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 57948, 0, 3,
                                                                       56958, 24046, 56988, 7792,
                                                                       7810, 24748, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58008, 0, 3,
                                                                       56988, 24064, 57018, 7810,
                                                                       7828, 24784, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58068, 0, 3,
                                                                       57048, 24136, 57078, 7864,
                                                                       7882, 24892, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58128, 0, 3,
                                                                       57078, 24154, 57108, 7882,
                                                                       7900, 24928, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58188, 0, 3,
                                                                       57108, 24172, 57138, 7900,
                                                                       7918, 24964, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58248, 0, 3,
                                                                       57138, 24190, 57168, 7918,
                                                                       7936, 25000, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58308, 0, 3,
                                                                       57168, 24208, 57198, 7936,
                                                                       7954, 25036, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58368, 0, 3,
                                                                       57198, 24226, 57228, 7954,
                                                                       7972, 25072, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58428, 0, 3,
                                                                       57228, 24244, 57258, 7972,
                                                                       7990, 25108, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58488, 0, 3,
                                                                       57258, 24262, 57288, 7990,
                                                                       8008, 25144, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58548, 0, 3,
                                                                       57288, 24280, 57318, 8008,
                                                                       8026, 25180, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58608, 0, 3,
                                                                       57318, 24298, 57348, 8026,
                                                                       8044, 25216, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58668, 0, 3,
                                                                       57348, 24316, 57378, 8044,
                                                                       8062, 25252, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 58728, 0, 3,
                                                                       57408, 24424, 57468, 8098,
                                                                       8128, 25408, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 58828, 0, 3,
                                                                       57468, 24460, 57528, 8128,
                                                                       8158, 25468, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 58928, 0, 3,
                                                                       57528, 24496, 57588, 8158,
                                                                       8188, 25528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 59028, 0, 3,
                                                                       57588, 24532, 57648, 8188,
                                                                       8218, 25588, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 59128, 0, 3,
                                                                       57648, 24568, 57708, 8218,
                                                                       8248, 25648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 59228, 0, 3,
                                                                       57708, 24604, 57768, 8248,
                                                                       8278, 25708, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 59328, 0, 3,
                                                                       57768, 24640, 57828, 8278,
                                                                       8308, 25768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 59428, 0, 3,
                                                                       57828, 24676, 57888, 8308,
                                                                       8338, 25828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 59528, 0, 3,
                                                                       57888, 24712, 57948, 8338,
                                                                       8368, 25888, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 59628, 0, 3,
                                                                       57948, 24748, 58008, 8368,
                                                                       8398, 25948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 59728, 0, 3,
                                                                       58068, 24892, 58128, 8458,
                                                                       8488, 26128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 59828, 0, 3,
                                                                       58128, 24928, 58188, 8488,
                                                                       8518, 26188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 59928, 0, 3,
                                                                       58188, 24964, 58248, 8518,
                                                                       8548, 26248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60028, 0, 3,
                                                                       58248, 25000, 58308, 8548,
                                                                       8578, 26308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60128, 0, 3,
                                                                       58308, 25036, 58368, 8578,
                                                                       8608, 26368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60228, 0, 3,
                                                                       58368, 25072, 58428, 8608,
                                                                       8638, 26428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60328, 0, 3,
                                                                       58428, 25108, 58488, 8638,
                                                                       8668, 26488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60428, 0, 3,
                                                                       58488, 25144, 58548, 8668,
                                                                       8698, 26548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60528, 0, 3,
                                                                       58548, 25180, 58608, 8698,
                                                                       8728, 26608, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60628, 0, 3,
                                                                       58608, 25216, 58668, 8728,
                                                                       8758, 26668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 60728, 0, 3,
                                                                       58728, 25408, 58828, 8818,
                                                                       8863, 26908, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 60878, 0, 3,
                                                                       58828, 25468, 58928, 8863,
                                                                       8908, 26998, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 61028, 0, 3,
                                                                       58928, 25528, 59028, 8908,
                                                                       8953, 27088, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 61178, 0, 3,
                                                                       59028, 25588, 59128, 8953,
                                                                       8998, 27178, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 61328, 0, 3,
                                                                       59128, 25648, 59228, 8998,
                                                                       9043, 27268, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 61478, 0, 3,
                                                                       59228, 25708, 59328, 9043,
                                                                       9088, 27358, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 61628, 0, 3,
                                                                       59328, 25768, 59428, 9088,
                                                                       9133, 27448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 61778, 0, 3,
                                                                       59428, 25828, 59528, 9133,
                                                                       9178, 27538, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 61928, 0, 3,
                                                                       59528, 25888, 59628, 9178,
                                                                       9223, 27628, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 62078, 0, 3,
                                                                       59728, 26128, 59828, 9313,
                                                                       9358, 27898, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 62228, 0, 3,
                                                                       59828, 26188, 59928, 9358,
                                                                       9403, 27988, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 62378, 0, 3,
                                                                       59928, 26248, 60028, 9403,
                                                                       9448, 28078, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 62528, 0, 3,
                                                                       60028, 26308, 60128, 9448,
                                                                       9493, 28168, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 62678, 0, 3,
                                                                       60128, 26368, 60228, 9493,
                                                                       9538, 28258, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 62828, 0, 3,
                                                                       60228, 26428, 60328, 9538,
                                                                       9583, 28348, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 62978, 0, 3,
                                                                       60328, 26488, 60428, 9583,
                                                                       9628, 28438, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 63128, 0, 3,
                                                                       60428, 26548, 60528, 9628,
                                                                       9673, 28528, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 63278, 0, 3,
                                                                       60528, 26608, 60628, 9673,
                                                                       9718, 28618, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 63428, 0, 3,
                                                                       60728, 26908, 60878, 9808,
                                                                       9871, 28960, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 63638, 0, 3,
                                                                       60878, 26998, 61028, 9871,
                                                                       9934, 29086, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 63848, 0, 3,
                                                                       61028, 27088, 61178, 9934,
                                                                       9997, 29212, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 64058, 0, 3,
                                                                       61178, 27178, 61328, 9997,
                                                                       10060, 29338, ncols,
                                                                       gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 64268, 0, 3,
                                                                       61328, 27268, 61478,
                                                                       10060, 10123, 29464,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 64478, 0, 3,
                                                                       61478, 27358, 61628,
                                                                       10123, 10186, 29590,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 64688, 0, 3,
                                                                       61628, 27448, 61778,
                                                                       10186, 10249, 29716,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 64898, 0, 3,
                                                                       61778, 27538, 61928,
                                                                       10249, 10312, 29842,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 65108, 0, 3,
                                                                       62078, 27898, 62228,
                                                                       10438, 10501, 30220,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 65318, 0, 3,
                                                                       62228, 27988, 62378,
                                                                       10501, 10564, 30346,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 65528, 0, 3,
                                                                       62378, 28078, 62528,
                                                                       10564, 10627, 30472,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 65738, 0, 3,
                                                                       62528, 28168, 62678,
                                                                       10627, 10690, 30598,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 65948, 0, 3,
                                                                       62678, 28258, 62828,
                                                                       10690, 10753, 30724,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 66158, 0, 3,
                                                                       62828, 28348, 62978,
                                                                       10753, 10816, 30850,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 66368, 0, 3,
                                                                       62978, 28438, 63128,
                                                                       10816, 10879, 30976,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 66578, 0, 3,
                                                                       63128, 28528, 63278,
                                                                       10879, 10942, 31102,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 66788, 0, 3,
                                                                       63428, 28960, 63638,
                                                                       11068, 11152, 31564,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 67068, 0, 3,
                                                                       63638, 29086, 63848,
                                                                       11152, 11236, 31732,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 67348, 0, 3,
                                                                       63848, 29212, 64058,
                                                                       11236, 11320, 31900,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 67628, 0, 3,
                                                                       64058, 29338, 64268,
                                                                       11320, 11404, 32068,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 67908, 0, 3,
                                                                       64268, 29464, 64478,
                                                                       11404, 11488, 32236,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 68188, 0, 3,
                                                                       64478, 29590, 64688,
                                                                       11488, 11572, 32404,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 68468, 0, 3,
                                                                       64688, 29716, 64898,
                                                                       11572, 11656, 32572,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 68748, 0, 3,
                                                                       65108, 30220, 65318,
                                                                       11824, 11908, 33076,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 69028, 0, 3,
                                                                       65318, 30346, 65528,
                                                                       11908, 11992, 33244,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 69308, 0, 3,
                                                                       65528, 30472, 65738,
                                                                       11992, 12076, 33412,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 69588, 0, 3,
                                                                       65738, 30598, 65948,
                                                                       12076, 12160, 33580,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 69868, 0, 3,
                                                                       65948, 30724, 66158,
                                                                       12160, 12244, 33748,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 70148, 0, 3,
                                                                       66158, 30850, 66368,
                                                                       12244, 12328, 33916,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 70428, 0, 3,
                                                                       66368, 30976, 66578,
                                                                       12328, 12412, 34084,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 70708, 0, 3,
                                                                       66788, 31564, 67068,
                                                                       12580, 12688, 34684,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 71068, 0, 3,
                                                                       67068, 31732, 67348,
                                                                       12688, 12796, 34900,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 71428, 0, 3,
                                                                       67348, 31900, 67628,
                                                                       12796, 12904, 35116,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 71788, 0, 3,
                                                                       67628, 32068, 67908,
                                                                       12904, 13012, 35332,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 72148, 0, 3,
                                                                       67908, 32236, 68188,
                                                                       13012, 13120, 35548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 72508, 0, 3,
                                                                       68188, 32404, 68468,
                                                                       13120, 13228, 35764,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 72868, 0, 3,
                                                                       68748, 33076, 69028,
                                                                       13444, 13552, 36412,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 73228, 0, 3,
                                                                       69028, 33244, 69308,
                                                                       13552, 13660, 36628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 73588, 0, 3,
                                                                       69308, 33412, 69588,
                                                                       13660, 13768, 36844,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 73948, 0, 3,
                                                                       69588, 33580, 69868,
                                                                       13768, 13876, 37060,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 74308, 0, 3,
                                                                       69868, 33748, 70148,
                                                                       13876, 13984, 37276,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 74668, 0, 3,
                                                                       70148, 33916, 70428,
                                                                       13984, 14092, 37492,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 75028, 0, 3,
                                                                       70708, 34684, 71068,
                                                                       14308, 14443, 38248,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 75478, 0, 3,
                                                                       71068, 34900, 71428,
                                                                       14443, 14578, 38518,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 75928, 0, 3,
                                                                       71428, 35116, 71788,
                                                                       14578, 14713, 38788,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 76378, 0, 3,
                                                                       71788, 35332, 72148,
                                                                       14713, 14848, 39058,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 76828, 0, 3,
                                                                       72148, 35548, 72508,
                                                                       14848, 14983, 39328,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 77278, 0, 3,
                                                                       72868, 36412, 73228,
                                                                       15253, 15388, 40138,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 77728, 0, 3,
                                                                       73228, 36628, 73588,
                                                                       15388, 15523, 40408,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 78178, 0, 3,
                                                                       73588, 36844, 73948,
                                                                       15523, 15658, 40678,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 78628, 0, 3,
                                                                       73948, 37060, 74308,
                                                                       15658, 15793, 40948,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 79078, 0, 3,
                                                                       74308, 37276, 74668,
                                                                       15793, 15928, 41218,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 79528, 0, 3,
                                                                       75028, 38248, 75478,
                                                                       16198, 16363, 42148,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 80078, 0, 3,
                                                                       75478, 38518, 75928,
                                                                       16363, 16528, 42478,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 80628, 0, 3,
                                                                       75928, 38788, 76378,
                                                                       16528, 16693, 42808,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 81178, 0, 3,
                                                                       76378, 39058, 76828,
                                                                       16693, 16858, 43138,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 81728, 0, 3,
                                                                       77278, 40138, 77728,
                                                                       17188, 17353, 44128,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 82278, 0, 3,
                                                                       77728, 40408, 78178,
                                                                       17353, 17518, 44458,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 82828, 0, 3,
                                                                       78178, 40678, 78628,
                                                                       17518, 17683, 44788,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 83378, 0, 3,
                                                                       78628, 40948, 79078,
                                                                       17683, 17848, 45118,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 83928, 0, 3,
                                                                       79528, 42148, 80078,
                                                                       18178, 18376, 46240,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 84588, 0, 3,
                                                                       80078, 42478, 80628,
                                                                       18376, 18574, 46636,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 85248, 0, 3,
                                                                       80628, 42808, 81178,
                                                                       18574, 18772, 47032,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 85908, 0, 3,
                                                                       81728, 44128, 82278,
                                                                       19168, 19366, 48220,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 86568, 0, 3,
                                                                       82278, 44458, 82828,
                                                                       19366, 19564, 48616,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 87228, 0, 3,
                                                                       82828, 44788, 83378,
                                                                       19564, 19762, 49012,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 87888, 0, 3,
                                                                       83928, 46240, 84588,
                                                                       20158, 20392, 50344,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 88668, 0, 3,
                                                                       84588, 46636, 85248,
                                                                       20392, 20626, 50812,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 89448, 0, 3,
                                                                       85908, 48220, 86568,
                                                                       21094, 21328, 52216,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 90228, 0, 3,
                                                                       86568, 48616, 87228,
                                                                       21328, 21562, 52684,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 91008, 0, 3,
                                                                       87888, 50344, 88668,
                                                                       22030, 22303, 54244,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 91918, 0, 3,
                                                                       89448, 52216, 90228,
                                                                       22849, 23122, 55882,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92828, 3, 23668,
                                                                       23674, 56428, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92843, 3, 23674,
                                                                       23680, 56438, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92858, 3, 23680,
                                                                       23686, 56448, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92873, 3, 23686,
                                                                       23692, 56458, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92888, 3, 23692,
                                                                       23698, 56468, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92903, 3, 23698,
                                                                       23704, 56478, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92918, 3, 23704,
                                                                       23710, 56488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92933, 3, 23710,
                                                                       23716, 56498, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92948, 3, 23716,
                                                                       23722, 56508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92963, 3, 23722,
                                                                       23728, 56518, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92978, 3, 23728,
                                                                       23734, 56528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 92993, 3, 23734,
                                                                       23740, 56538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93008, 3, 23740,
                                                                       23746, 56548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93023, 3, 23758,
                                                                       23764, 56558, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93038, 3, 23764,
                                                                       23770, 56568, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93053, 3, 23770,
                                                                       23776, 56578, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93068, 3, 23776,
                                                                       23782, 56588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93083, 3, 23782,
                                                                       23788, 56598, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93098, 3, 23788,
                                                                       23794, 56608, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93113, 3, 23794,
                                                                       23800, 56618, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93128, 3, 23800,
                                                                       23806, 56628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93143, 3, 23806,
                                                                       23812, 56638, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93158, 3, 23812,
                                                                       23818, 56648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93173, 3, 23818,
                                                                       23824, 56658, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93188, 3, 23824,
                                                                       23830, 56668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 93203, 3, 23830,
                                                                       23836, 56678, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93218, 0, 3,
                                                                       92828, 56428, 92843,
                                                                       23848, 23866, 56688,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93263, 0, 3,
                                                                       92843, 56438, 92858,
                                                                       23866, 23884, 56718,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93308, 0, 3,
                                                                       92858, 56448, 92873,
                                                                       23884, 23902, 56748,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93353, 0, 3,
                                                                       92873, 56458, 92888,
                                                                       23902, 23920, 56778,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93398, 0, 3,
                                                                       92888, 56468, 92903,
                                                                       23920, 23938, 56808,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93443, 0, 3,
                                                                       92903, 56478, 92918,
                                                                       23938, 23956, 56838,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93488, 0, 3,
                                                                       92918, 56488, 92933,
                                                                       23956, 23974, 56868,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93533, 0, 3,
                                                                       92933, 56498, 92948,
                                                                       23974, 23992, 56898,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93578, 0, 3,
                                                                       92948, 56508, 92963,
                                                                       23992, 24010, 56928,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93623, 0, 3,
                                                                       92963, 56518, 92978,
                                                                       24010, 24028, 56958,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93668, 0, 3,
                                                                       92978, 56528, 92993,
                                                                       24028, 24046, 56988,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93713, 0, 3,
                                                                       92993, 56538, 93008,
                                                                       24046, 24064, 57018,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93758, 0, 3,
                                                                       93023, 56558, 93038,
                                                                       24100, 24118, 57048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93803, 0, 3,
                                                                       93038, 56568, 93053,
                                                                       24118, 24136, 57078,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93848, 0, 3,
                                                                       93053, 56578, 93068,
                                                                       24136, 24154, 57108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93893, 0, 3,
                                                                       93068, 56588, 93083,
                                                                       24154, 24172, 57138,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93938, 0, 3,
                                                                       93083, 56598, 93098,
                                                                       24172, 24190, 57168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 93983, 0, 3,
                                                                       93098, 56608, 93113,
                                                                       24190, 24208, 57198,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 94028, 0, 3,
                                                                       93113, 56618, 93128,
                                                                       24208, 24226, 57228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 94073, 0, 3,
                                                                       93128, 56628, 93143,
                                                                       24226, 24244, 57258,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 94118, 0, 3,
                                                                       93143, 56638, 93158,
                                                                       24244, 24262, 57288,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 94163, 0, 3,
                                                                       93158, 56648, 93173,
                                                                       24262, 24280, 57318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 94208, 0, 3,
                                                                       93173, 56658, 93188,
                                                                       24280, 24298, 57348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 94253, 0, 3,
                                                                       93188, 56668, 93203,
                                                                       24298, 24316, 57378,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 94298, 0, 3,
                                                                       93218, 56688, 93263,
                                                                       24352, 24388, 57408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 94388, 0, 3,
                                                                       93263, 56718, 93308,
                                                                       24388, 24424, 57468,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 94478, 0, 3,
                                                                       93308, 56748, 93353,
                                                                       24424, 24460, 57528,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 94568, 0, 3,
                                                                       93353, 56778, 93398,
                                                                       24460, 24496, 57588,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 94658, 0, 3,
                                                                       93398, 56808, 93443,
                                                                       24496, 24532, 57648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 94748, 0, 3,
                                                                       93443, 56838, 93488,
                                                                       24532, 24568, 57708,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 94838, 0, 3,
                                                                       93488, 56868, 93533,
                                                                       24568, 24604, 57768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 94928, 0, 3,
                                                                       93533, 56898, 93578,
                                                                       24604, 24640, 57828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 95018, 0, 3,
                                                                       93578, 56928, 93623,
                                                                       24640, 24676, 57888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 95108, 0, 3,
                                                                       93623, 56958, 93668,
                                                                       24676, 24712, 57948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 95198, 0, 3,
                                                                       93668, 56988, 93713,
                                                                       24712, 24748, 58008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 95288, 0, 3,
                                                                       93758, 57048, 93803,
                                                                       24820, 24856, 58068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 95378, 0, 3,
                                                                       93803, 57078, 93848,
                                                                       24856, 24892, 58128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 95468, 0, 3,
                                                                       93848, 57108, 93893,
                                                                       24892, 24928, 58188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 95558, 0, 3,
                                                                       93893, 57138, 93938,
                                                                       24928, 24964, 58248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 95648, 0, 3,
                                                                       93938, 57168, 93983,
                                                                       24964, 25000, 58308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 95738, 0, 3,
                                                                       93983, 57198, 94028,
                                                                       25000, 25036, 58368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 95828, 0, 3,
                                                                       94028, 57228, 94073,
                                                                       25036, 25072, 58428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 95918, 0, 3,
                                                                       94073, 57258, 94118,
                                                                       25072, 25108, 58488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 96008, 0, 3,
                                                                       94118, 57288, 94163,
                                                                       25108, 25144, 58548,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 96098, 0, 3,
                                                                       94163, 57318, 94208,
                                                                       25144, 25180, 58608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 96188, 0, 3,
                                                                       94208, 57348, 94253,
                                                                       25180, 25216, 58668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 96278, 0, 3,
                                                                       94298, 57408, 94388,
                                                                       25288, 25348, 58728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 96428, 0, 3,
                                                                       94388, 57468, 94478,
                                                                       25348, 25408, 58828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 96578, 0, 3,
                                                                       94478, 57528, 94568,
                                                                       25408, 25468, 58928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 96728, 0, 3,
                                                                       94568, 57588, 94658,
                                                                       25468, 25528, 59028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 96878, 0, 3,
                                                                       94658, 57648, 94748,
                                                                       25528, 25588, 59128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 97028, 0, 3,
                                                                       94748, 57708, 94838,
                                                                       25588, 25648, 59228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 97178, 0, 3,
                                                                       94838, 57768, 94928,
                                                                       25648, 25708, 59328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 97328, 0, 3,
                                                                       94928, 57828, 95018,
                                                                       25708, 25768, 59428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 97478, 0, 3,
                                                                       95018, 57888, 95108,
                                                                       25768, 25828, 59528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 97628, 0, 3,
                                                                       95108, 57948, 95198,
                                                                       25828, 25888, 59628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 97778, 0, 3,
                                                                       95288, 58068, 95378,
                                                                       26008, 26068, 59728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 97928, 0, 3,
                                                                       95378, 58128, 95468,
                                                                       26068, 26128, 59828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 98078, 0, 3,
                                                                       95468, 58188, 95558,
                                                                       26128, 26188, 59928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 98228, 0, 3,
                                                                       95558, 58248, 95648,
                                                                       26188, 26248, 60028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 98378, 0, 3,
                                                                       95648, 58308, 95738,
                                                                       26248, 26308, 60128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 98528, 0, 3,
                                                                       95738, 58368, 95828,
                                                                       26308, 26368, 60228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 98678, 0, 3,
                                                                       95828, 58428, 95918,
                                                                       26368, 26428, 60328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 98828, 0, 3,
                                                                       95918, 58488, 96008,
                                                                       26428, 26488, 60428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 98978, 0, 3,
                                                                       96008, 58548, 96098,
                                                                       26488, 26548, 60528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 99128, 0, 3,
                                                                       96098, 58608, 96188,
                                                                       26548, 26608, 60628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 99278, 0, 3,
                                                                       96278, 58728, 96428,
                                                                       26728, 26818, 60728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 99503, 0, 3,
                                                                       96428, 58828, 96578,
                                                                       26818, 26908, 60878,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 99728, 0, 3,
                                                                       96578, 58928, 96728,
                                                                       26908, 26998, 61028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 99953, 0, 3,
                                                                       96728, 59028, 96878,
                                                                       26998, 27088, 61178,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 100178, 0, 3,
                                                                       96878, 59128, 97028,
                                                                       27088, 27178, 61328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 100403, 0, 3,
                                                                       97028, 59228, 97178,
                                                                       27178, 27268, 61478,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 100628, 0, 3,
                                                                       97178, 59328, 97328,
                                                                       27268, 27358, 61628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 100853, 0, 3,
                                                                       97328, 59428, 97478,
                                                                       27358, 27448, 61778,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 101078, 0, 3,
                                                                       97478, 59528, 97628,
                                                                       27448, 27538, 61928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 101303, 0, 3,
                                                                       97778, 59728, 97928,
                                                                       27718, 27808, 62078,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 101528, 0, 3,
                                                                       97928, 59828, 98078,
                                                                       27808, 27898, 62228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 101753, 0, 3,
                                                                       98078, 59928, 98228,
                                                                       27898, 27988, 62378,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 101978, 0, 3,
                                                                       98228, 60028, 98378,
                                                                       27988, 28078, 62528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 102203, 0, 3,
                                                                       98378, 60128, 98528,
                                                                       28078, 28168, 62678,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 102428, 0, 3,
                                                                       98528, 60228, 98678,
                                                                       28168, 28258, 62828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 102653, 0, 3,
                                                                       98678, 60328, 98828,
                                                                       28258, 28348, 62978,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 102878, 0, 3,
                                                                       98828, 60428, 98978,
                                                                       28348, 28438, 63128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 103103, 0, 3,
                                                                       98978, 60528, 99128,
                                                                       28438, 28528, 63278,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 103328, 0, 3,
                                                                       99278, 60728, 99503,
                                                                       28708, 28834, 63428,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 103643, 0, 3,
                                                                       99503, 60878, 99728,
                                                                       28834, 28960, 63638,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 103958, 0, 3,
                                                                       99728, 61028, 99953,
                                                                       28960, 29086, 63848,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 104273, 0, 3,
                                                                       99953, 61178, 100178,
                                                                       29086, 29212, 64058,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 104588, 0, 3,
                                                                       100178, 61328, 100403,
                                                                       29212, 29338, 64268,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 104903, 0, 3,
                                                                       100403, 61478, 100628,
                                                                       29338, 29464, 64478,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 105218, 0, 3,
                                                                       100628, 61628, 100853,
                                                                       29464, 29590, 64688,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 105533, 0, 3,
                                                                       100853, 61778, 101078,
                                                                       29590, 29716, 64898,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 105848, 0, 3,
                                                                       101303, 62078, 101528,
                                                                       29968, 30094, 65108,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 106163, 0, 3,
                                                                       101528, 62228, 101753,
                                                                       30094, 30220, 65318,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 106478, 0, 3,
                                                                       101753, 62378, 101978,
                                                                       30220, 30346, 65528,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 106793, 0, 3,
                                                                       101978, 62528, 102203,
                                                                       30346, 30472, 65738,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 107108, 0, 3,
                                                                       102203, 62678, 102428,
                                                                       30472, 30598, 65948,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 107423, 0, 3,
                                                                       102428, 62828, 102653,
                                                                       30598, 30724, 66158,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 107738, 0, 3,
                                                                       102653, 62978, 102878,
                                                                       30724, 30850, 66368,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 108053, 0, 3,
                                                                       102878, 63128, 103103,
                                                                       30850, 30976, 66578,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 108368, 0, 3,
                                                                       103328, 63428, 103643,
                                                                       31228, 31396, 66788,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 108788, 0, 3,
                                                                       103643, 63638, 103958,
                                                                       31396, 31564, 67068,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 109208, 0, 3,
                                                                       103958, 63848, 104273,
                                                                       31564, 31732, 67348,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 109628, 0, 3,
                                                                       104273, 64058, 104588,
                                                                       31732, 31900, 67628,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 110048, 0, 3,
                                                                       104588, 64268, 104903,
                                                                       31900, 32068, 67908,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 110468, 0, 3,
                                                                       104903, 64478, 105218,
                                                                       32068, 32236, 68188,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 110888, 0, 3,
                                                                       105218, 64688, 105533,
                                                                       32236, 32404, 68468,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 111308, 0, 3,
                                                                       105848, 65108, 106163,
                                                                       32740, 32908, 68748,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 111728, 0, 3,
                                                                       106163, 65318, 106478,
                                                                       32908, 33076, 69028,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 112148, 0, 3,
                                                                       106478, 65528, 106793,
                                                                       33076, 33244, 69308,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 112568, 0, 3,
                                                                       106793, 65738, 107108,
                                                                       33244, 33412, 69588,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 112988, 0, 3,
                                                                       107108, 65948, 107423,
                                                                       33412, 33580, 69868,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 113408, 0, 3,
                                                                       107423, 66158, 107738,
                                                                       33580, 33748, 70148,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 113828, 0, 3,
                                                                       107738, 66368, 108053,
                                                                       33748, 33916, 70428,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 114248, 0, 3,
                                                                       108368, 66788, 108788,
                                                                       34252, 34468, 70708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 114788, 0, 3,
                                                                       108788, 67068, 109208,
                                                                       34468, 34684, 71068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 115328, 0, 3,
                                                                       109208, 67348, 109628,
                                                                       34684, 34900, 71428,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 115868, 0, 3,
                                                                       109628, 67628, 110048,
                                                                       34900, 35116, 71788,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 116408, 0, 3,
                                                                       110048, 67908, 110468,
                                                                       35116, 35332, 72148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 116948, 0, 3,
                                                                       110468, 68188, 110888,
                                                                       35332, 35548, 72508,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 117488, 0, 3,
                                                                       111308, 68748, 111728,
                                                                       35980, 36196, 72868,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 118028, 0, 3,
                                                                       111728, 69028, 112148,
                                                                       36196, 36412, 73228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 118568, 0, 3,
                                                                       112148, 69308, 112568,
                                                                       36412, 36628, 73588,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 119108, 0, 3,
                                                                       112568, 69588, 112988,
                                                                       36628, 36844, 73948,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 119648, 0, 3,
                                                                       112988, 69868, 113408,
                                                                       36844, 37060, 74308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 120188, 0, 3,
                                                                       113408, 70148, 113828,
                                                                       37060, 37276, 74668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 120728, 0, 3,
                                                                       114248, 70708, 114788,
                                                                       37708, 37978, 75028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 121403, 0, 3,
                                                                       114788, 71068, 115328,
                                                                       37978, 38248, 75478,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 122078, 0, 3,
                                                                       115328, 71428, 115868,
                                                                       38248, 38518, 75928,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 122753, 0, 3,
                                                                       115868, 71788, 116408,
                                                                       38518, 38788, 76378,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 123428, 0, 3,
                                                                       116408, 72148, 116948,
                                                                       38788, 39058, 76828,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 124103, 0, 3,
                                                                       117488, 72868, 118028,
                                                                       39598, 39868, 77278,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 124778, 0, 3,
                                                                       118028, 73228, 118568,
                                                                       39868, 40138, 77728,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 125453, 0, 3,
                                                                       118568, 73588, 119108,
                                                                       40138, 40408, 78178,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 126128, 0, 3,
                                                                       119108, 73948, 119648,
                                                                       40408, 40678, 78628,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 126803, 0, 3,
                                                                       119648, 74308, 120188,
                                                                       40678, 40948, 79078,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 127478, 0, 3,
                                                                       120728, 75028, 121403,
                                                                       41488, 41818, 79528,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 128303, 0, 3,
                                                                       121403, 75478, 122078,
                                                                       41818, 42148, 80078,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 129128, 0, 3,
                                                                       122078, 75928, 122753,
                                                                       42148, 42478, 80628,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 129953, 0, 3,
                                                                       122753, 76378, 123428,
                                                                       42478, 42808, 81178,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 130778, 0, 3,
                                                                       124103, 77278, 124778,
                                                                       43468, 43798, 81728,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 131603, 0, 3,
                                                                       124778, 77728, 125453,
                                                                       43798, 44128, 82278,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 132428, 0, 3,
                                                                       125453, 78178, 126128,
                                                                       44128, 44458, 82828,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 133253, 0, 3,
                                                                       126128, 78628, 126803,
                                                                       44458, 44788, 83378,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 134078, 0, 3,
                                                                       127478, 79528, 128303,
                                                                       45448, 45844, 83928,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 135068, 0, 3,
                                                                       128303, 80078, 129128,
                                                                       45844, 46240, 84588,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 136058, 0, 3,
                                                                       129128, 80628, 129953,
                                                                       46240, 46636, 85248,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 137048, 0, 3,
                                                                       130778, 81728, 131603,
                                                                       47428, 47824, 85908,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 138038, 0, 3,
                                                                       131603, 82278, 132428,
                                                                       47824, 48220, 86568,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 139028, 0, 3,
                                                                       132428, 82828, 133253,
                                                                       48220, 48616, 87228,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 140018, 0, 3,
                                                                       134078, 83928, 135068,
                                                                       49408, 49876, 87888,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 141188, 0, 3,
                                                                       135068, 84588, 136058,
                                                                       49876, 50344, 88668,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 142358, 0, 3,
                                                                       137048, 85908, 138038,
                                                                       51280, 51748, 89448,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 143528, 0, 3,
                                                                       138038, 86568, 139028,
                                                                       51748, 52216, 90228,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 144698, 0, 3,
                                                                       140018, 87888, 141188,
                                                                       53152, 53698, 91008,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 146063, 0, 3,
                                                                       142358, 89448, 143528,
                                                                       54790, 55336, 91918,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 147428, 108368, 420, ncols);

                    simdfunc::contract_primitives(buffer, 148100, 111308, 420, ncols);

                    simdfunc::contract_primitives(buffer, 148772, 114248, 540, ncols);

                    simdfunc::contract_primitives(buffer, 149636, 117488, 540, ncols);

                    simdfunc::contract_primitives(buffer, 150500, 120728, 675, ncols);

                    simdfunc::contract_primitives(buffer, 151580, 124103, 675, ncols);

                    simdfunc::contract_primitives(buffer, 152660, 127478, 825, ncols);

                    simdfunc::contract_primitives(buffer, 153980, 130778, 825, ncols);

                    simdfunc::contract_primitives(buffer, 155300, 134078, 990, ncols);

                    simdfunc::contract_primitives(buffer, 156884, 137048, 990, ncols);

                    simdfunc::contract_primitives(buffer, 158468, 140018, 1170, ncols);

                    simdfunc::contract_primitives(buffer, 160340, 142358, 1170, ncols);

                    simdfunc::contract_primitives(buffer, 162212, 144698, 1365, ncols);

                    simdfunc::contract_primitives(buffer, 164396, 146063, 1365, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 147848, 147428, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 148520, 148100, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 149312, 148772, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 150176, 149636, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 151175, 150500, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 152255, 151580, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 153485, 152660, 55, 1, nmax);

        simdtrf::transform_g_inner(buffer, 154805, 153980, 55, 1, nmax);

        simdtrf::transform_g_inner(buffer, 156290, 155300, 66, 1, nmax);

        simdtrf::transform_g_inner(buffer, 157874, 156884, 66, 1, nmax);

        simdtrf::transform_g_inner(buffer, 159638, 158468, 78, 1, nmax);

        simdtrf::transform_g_inner(buffer, 161510, 160340, 78, 1, nmax);

        simdtrf::transform_g_inner(buffer, 163577, 162212, 91, 1, nmax);

        simdtrf::transform_g_inner(buffer, 165761, 164396, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 166580, 147848, 149312, 9,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 167336, 148520, 150176, 9,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 168092, 149312, 151175, 9,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 169064, 150176, 152255, 9,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 170036, 151175, 153485, 9,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 171251, 152255, 154805, 9,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 172466, 153485, 156290, 9,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 173951, 154805, 157874, 9,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 175436, 156290, 159638, 9,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 177218, 157874, 161510, 9,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 179000, 159638, 163577, 9,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 181106, 161510, 165761, 9,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 183212, 166580, 168092, 9,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 184724, 167336, 169064, 9,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 186236, 168092, 170036, 9,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 188180, 169064, 171251, 9,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 190124, 170036, 172466, 9,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 192554, 171251, 173951, 9,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 194984, 172466, 175436, 9,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 197954, 173951, 177218, 9,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 200924, 175436, 179000, 9,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 204488, 177218, 181106, 9,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 208052, 183212, 186236, 9,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 210572, 184724, 188180, 9,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 213092, 186236, 190124, 9,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 216332, 188180, 192554, 9,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 219572, 190124, 194984, 9,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 223622, 192554, 197954, 9,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 227672, 194984, 200924, 9,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 232622, 197954, 204488, 9,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 237572, 208052, 213092, 9,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 241352, 210572, 216332, 9,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 245132, 213092, 219572, 9,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 249992, 216332, 223622, 9,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 254852, 219572, 227672, 9,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 260927, 223622, 232622, 9,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 267002, 237572, 245132, 9,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 272294, 241352, 249992, 9,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 277586, 245132, 254852, 9,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 284390, 249992, 260927, 9,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 291194, 267002, 277586, 9,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 298250, 272294, 284390, 9,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 305306, 298250, 28, 9, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 305306, 117, nmax);

        simdtrf::transform_i_inner(buffer, 305306, 291194, 28, 9, nmax);

        simdtrf::transform_i_outer(values + 1521 * nvalues + n * npairs, nvalues, buffer, 305306,
                                   117, nmax);
    }

    for (size_t m = 0; m < 3042; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
