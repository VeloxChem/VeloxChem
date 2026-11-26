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


#include "SimdThreeCenterElectronRepulsionRsRecHHI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
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
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_hhi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hhi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 331771, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3146 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 331771, 220800, 19724, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5442, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5445, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5448, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5451, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5454, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5457, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5460, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5463, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5466, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5469, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5472, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5475, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5478, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5481, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5484, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5487, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5490, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5493, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5496, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5499, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5502, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5505, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5508, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5511, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5514, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5517, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5520, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5523, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5526, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5529, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5532, 3, 9, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5541, 3, 10, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5550, 3, 11, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5559, 3, 12, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5568, 3, 13, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5577, 3, 14, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5586, 3, 15, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5595, 3, 16, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5604, 3, 17, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5613, 3, 18, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5622, 3, 19, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5631, 3, 20, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5640, 3, 21, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5649, 3, 22, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5658, 3, 27, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5667, 3, 28, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5676, 3, 29, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5685, 3, 30, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5694, 3, 31, 108,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5703, 3, 32, 111,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5712, 3, 33, 114,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5721, 3, 34, 117,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5730, 3, 35, 120,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5739, 3, 36, 123,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5748, 3, 37, 126,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5757, 3, 38, 129,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5766, 3, 39, 132,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5775, 3, 40, 135,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5784, 3, 48, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5802, 3, 51, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5820, 3, 54, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5838, 3, 57, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5856, 3, 60, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5874, 3, 63, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5892, 3, 66, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5910, 3, 69, 192,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5928, 3, 72, 198,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5946, 3, 75, 204,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5964, 3, 78, 210,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5982, 3, 81, 216,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6000, 3, 84, 222,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6018, 3, 96, 240,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6036, 3, 99, 246,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6054, 3, 102, 252,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6072, 3, 105, 258,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6090, 3, 108, 264,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6108, 3, 111, 270,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6126, 3, 114, 276,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6144, 3, 117, 282,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6162, 3, 120, 288,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6180, 3, 123, 294,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6198, 3, 126, 300,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6216, 3, 129, 306,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6234, 3, 132, 312,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6252, 3, 150, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6282, 3, 156, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6312, 3, 162, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6342, 3, 168, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6372, 3, 174, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6402, 3, 180, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6432, 3, 186, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6462, 3, 192, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6492, 3, 198, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6522, 3, 204, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6552, 3, 210, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6582, 3, 216, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6612, 3, 240, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6642, 3, 246, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6672, 3, 252, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6702, 3, 258, 508,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6732, 3, 264, 518,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6762, 3, 270, 528,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6792, 3, 276, 538,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6822, 3, 282, 548,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6852, 3, 288, 558,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6882, 3, 294, 568,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6912, 3, 300, 578,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6942, 3, 306, 588,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6972, 3, 338, 628,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7017, 3, 348, 643,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7062, 3, 358, 658,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7107, 3, 368, 673,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7152, 3, 378, 688,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7197, 3, 388, 703,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7242, 3, 398, 718,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7287, 3, 408, 733,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7332, 3, 418, 748,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7377, 3, 428, 763,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7422, 3, 438, 778,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7467, 3, 478, 823,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7512, 3, 488, 838,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7557, 3, 498, 853,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7602, 3, 508, 868,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7647, 3, 518, 883,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7692, 3, 528, 898,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7737, 3, 538, 913,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7782, 3, 548, 928,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7827, 3, 558, 943,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7872, 3, 568, 958,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7917, 3, 578, 973,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7962, 3, 628,
                                                                       1030, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8025, 3, 643,
                                                                       1051, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8088, 3, 658,
                                                                       1072, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8151, 3, 673,
                                                                       1093, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8214, 3, 688,
                                                                       1114, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8277, 3, 703,
                                                                       1135, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8340, 3, 718,
                                                                       1156, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8403, 3, 733,
                                                                       1177, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8466, 3, 748,
                                                                       1198, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8529, 3, 763,
                                                                       1219, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8592, 3, 823,
                                                                       1282, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8655, 3, 838,
                                                                       1303, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8718, 3, 853,
                                                                       1324, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8781, 3, 868,
                                                                       1345, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8844, 3, 883,
                                                                       1366, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8907, 3, 898,
                                                                       1387, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8970, 3, 913,
                                                                       1408, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9033, 3, 928,
                                                                       1429, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9096, 3, 943,
                                                                       1450, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9159, 3, 958,
                                                                       1471, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9222, 3, 1030,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9306, 3, 1051,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9390, 3, 1072,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9474, 3, 1093,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9558, 3, 1114,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9642, 3, 1135,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9726, 3, 1156,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9810, 3, 1177,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9894, 3, 1198,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9978, 3, 1282,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10062, 3, 1303,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10146, 3, 1324,
                                                                       1912, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10230, 3, 1345,
                                                                       1940, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10314, 3, 1366,
                                                                       1968, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10398, 3, 1387,
                                                                       1996, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10482, 3, 1408,
                                                                       2024, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10566, 3, 1429,
                                                                       2052, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10650, 3, 1450,
                                                                       2080, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10734, 3, 1548,
                                                                       2180, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10842, 3, 1576,
                                                                       2216, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10950, 3, 1604,
                                                                       2252, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11058, 3, 1632,
                                                                       2288, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11166, 3, 1660,
                                                                       2324, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11274, 3, 1688,
                                                                       2360, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11382, 3, 1716,
                                                                       2396, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11490, 3, 1744,
                                                                       2432, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11598, 3, 1856,
                                                                       2540, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11706, 3, 1884,
                                                                       2576, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11814, 3, 1912,
                                                                       2612, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11922, 3, 1940,
                                                                       2648, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12030, 3, 1968,
                                                                       2684, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12138, 3, 1996,
                                                                       2720, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12246, 3, 2024,
                                                                       2756, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12354, 3, 2052,
                                                                       2792, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12462, 3, 2180,
                                                                       2918, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12597, 3, 2216,
                                                                       2963, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12732, 3, 2252,
                                                                       3008, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12867, 3, 2288,
                                                                       3053, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13002, 3, 2324,
                                                                       3098, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13137, 3, 2360,
                                                                       3143, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13272, 3, 2396,
                                                                       3188, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13407, 3, 2540,
                                                                       3323, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13542, 3, 2576,
                                                                       3368, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13677, 3, 2612,
                                                                       3413, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13812, 3, 2648,
                                                                       3458, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13947, 3, 2684,
                                                                       3503, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14082, 3, 2720,
                                                                       3548, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14217, 3, 2756,
                                                                       3593, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14352, 3, 2918,
                                                                       3748, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14517, 3, 2963,
                                                                       3803, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14682, 3, 3008,
                                                                       3858, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14847, 3, 3053,
                                                                       3913, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15012, 3, 3098,
                                                                       3968, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15177, 3, 3143,
                                                                       4023, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15342, 3, 3323,
                                                                       4188, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15507, 3, 3368,
                                                                       4243, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15672, 3, 3413,
                                                                       4298, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15837, 3, 3458,
                                                                       4353, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16002, 3, 3503,
                                                                       4408, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16167, 3, 3548,
                                                                       4463, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16332, 3, 3748,
                                                                       4650, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16530, 3, 3803,
                                                                       4716, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16728, 3, 3858,
                                                                       4782, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16926, 3, 3913,
                                                                       4848, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 17124, 3, 3968,
                                                                       4914, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 17322, 3, 4188,
                                                                       5112, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 17520, 3, 4243,
                                                                       5178, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 17718, 3, 4298,
                                                                       5244, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 17916, 3, 4353,
                                                                       5310, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18114, 3, 4408,
                                                                       5376, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18312, 3, 7, 8,
                                                                       5442, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18318, 3, 8, 9,
                                                                       5445, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18324, 3, 9, 10,
                                                                       5448, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18330, 3, 10, 11,
                                                                       5451, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18336, 3, 11, 12,
                                                                       5454, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18342, 3, 12, 13,
                                                                       5457, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18348, 3, 13, 14,
                                                                       5460, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18354, 3, 14, 15,
                                                                       5463, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18360, 3, 15, 16,
                                                                       5466, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18366, 3, 16, 17,
                                                                       5469, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18372, 3, 17, 18,
                                                                       5472, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18378, 3, 18, 19,
                                                                       5475, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18384, 3, 19, 20,
                                                                       5478, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18390, 3, 20, 21,
                                                                       5481, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18396, 3, 21, 22,
                                                                       5484, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18402, 3, 25, 26,
                                                                       5487, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18408, 3, 26, 27,
                                                                       5490, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18414, 3, 27, 28,
                                                                       5493, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18420, 3, 28, 29,
                                                                       5496, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18426, 3, 29, 30,
                                                                       5499, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18432, 3, 30, 31,
                                                                       5502, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18438, 3, 31, 32,
                                                                       5505, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18444, 3, 32, 33,
                                                                       5508, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18450, 3, 33, 34,
                                                                       5511, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18456, 3, 34, 35,
                                                                       5514, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18462, 3, 35, 36,
                                                                       5517, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18468, 3, 36, 37,
                                                                       5520, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18474, 3, 37, 38,
                                                                       5523, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18480, 3, 38, 39,
                                                                       5526, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18486, 3, 39, 40,
                                                                       5529, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18492, 0, 3,
                                                                       18312, 5442, 18318, 5532,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18510, 0, 3,
                                                                       18318, 5445, 18324, 5541,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18528, 0, 3,
                                                                       18324, 5448, 18330, 5550,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18546, 0, 3,
                                                                       18330, 5451, 18336, 5559,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18564, 0, 3,
                                                                       18336, 5454, 18342, 5568,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18582, 0, 3,
                                                                       18342, 5457, 18348, 5577,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18600, 0, 3,
                                                                       18348, 5460, 18354, 5586,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18618, 0, 3,
                                                                       18354, 5463, 18360, 5595,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18636, 0, 3,
                                                                       18360, 5466, 18366, 5604,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18654, 0, 3,
                                                                       18366, 5469, 18372, 5613,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18672, 0, 3,
                                                                       18372, 5472, 18378, 5622,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18690, 0, 3,
                                                                       18378, 5475, 18384, 5631,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18708, 0, 3,
                                                                       18384, 5478, 18390, 5640,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18726, 0, 3,
                                                                       18390, 5481, 18396, 5649,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18744, 0, 3,
                                                                       18402, 5487, 18408, 5658,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18762, 0, 3,
                                                                       18408, 5490, 18414, 5667,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18780, 0, 3,
                                                                       18414, 5493, 18420, 5676,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18798, 0, 3,
                                                                       18420, 5496, 18426, 5685,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18816, 0, 3,
                                                                       18426, 5499, 18432, 5694,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18834, 0, 3,
                                                                       18432, 5502, 18438, 5703,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18852, 0, 3,
                                                                       18438, 5505, 18444, 5712,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18870, 0, 3,
                                                                       18444, 5508, 18450, 5721,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18888, 0, 3,
                                                                       18450, 5511, 18456, 5730,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18906, 0, 3,
                                                                       18456, 5514, 18462, 5739,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18924, 0, 3,
                                                                       18462, 5517, 18468, 5748,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18942, 0, 3,
                                                                       18468, 5520, 18474, 5757,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18960, 0, 3,
                                                                       18474, 5523, 18480, 5766,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18978, 0, 3,
                                                                       18480, 5526, 18486, 5775,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18996, 0, 3,
                                                                       18492, 5532, 18510, 138,
                                                                       144, 5784, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19032, 0, 3,
                                                                       18510, 5541, 18528, 144,
                                                                       150, 5802, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19068, 0, 3,
                                                                       18528, 5550, 18546, 150,
                                                                       156, 5820, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19104, 0, 3,
                                                                       18546, 5559, 18564, 156,
                                                                       162, 5838, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19140, 0, 3,
                                                                       18564, 5568, 18582, 162,
                                                                       168, 5856, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19176, 0, 3,
                                                                       18582, 5577, 18600, 168,
                                                                       174, 5874, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19212, 0, 3,
                                                                       18600, 5586, 18618, 174,
                                                                       180, 5892, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19248, 0, 3,
                                                                       18618, 5595, 18636, 180,
                                                                       186, 5910, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19284, 0, 3,
                                                                       18636, 5604, 18654, 186,
                                                                       192, 5928, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19320, 0, 3,
                                                                       18654, 5613, 18672, 192,
                                                                       198, 5946, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19356, 0, 3,
                                                                       18672, 5622, 18690, 198,
                                                                       204, 5964, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19392, 0, 3,
                                                                       18690, 5631, 18708, 204,
                                                                       210, 5982, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19428, 0, 3,
                                                                       18708, 5640, 18726, 210,
                                                                       216, 6000, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19464, 0, 3,
                                                                       18744, 5658, 18762, 228,
                                                                       234, 6018, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19500, 0, 3,
                                                                       18762, 5667, 18780, 234,
                                                                       240, 6036, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19536, 0, 3,
                                                                       18780, 5676, 18798, 240,
                                                                       246, 6054, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19572, 0, 3,
                                                                       18798, 5685, 18816, 246,
                                                                       252, 6072, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19608, 0, 3,
                                                                       18816, 5694, 18834, 252,
                                                                       258, 6090, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19644, 0, 3,
                                                                       18834, 5703, 18852, 258,
                                                                       264, 6108, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19680, 0, 3,
                                                                       18852, 5712, 18870, 264,
                                                                       270, 6126, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19716, 0, 3,
                                                                       18870, 5721, 18888, 270,
                                                                       276, 6144, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19752, 0, 3,
                                                                       18888, 5730, 18906, 276,
                                                                       282, 6162, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19788, 0, 3,
                                                                       18906, 5739, 18924, 282,
                                                                       288, 6180, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19824, 0, 3,
                                                                       18924, 5748, 18942, 288,
                                                                       294, 6198, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19860, 0, 3,
                                                                       18942, 5757, 18960, 294,
                                                                       300, 6216, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19896, 0, 3,
                                                                       18960, 5766, 18978, 300,
                                                                       306, 6234, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19932, 0, 3,
                                                                       18996, 5784, 19032, 318,
                                                                       328, 6252, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19992, 0, 3,
                                                                       19032, 5802, 19068, 328,
                                                                       338, 6282, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20052, 0, 3,
                                                                       19068, 5820, 19104, 338,
                                                                       348, 6312, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20112, 0, 3,
                                                                       19104, 5838, 19140, 348,
                                                                       358, 6342, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20172, 0, 3,
                                                                       19140, 5856, 19176, 358,
                                                                       368, 6372, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20232, 0, 3,
                                                                       19176, 5874, 19212, 368,
                                                                       378, 6402, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20292, 0, 3,
                                                                       19212, 5892, 19248, 378,
                                                                       388, 6432, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20352, 0, 3,
                                                                       19248, 5910, 19284, 388,
                                                                       398, 6462, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20412, 0, 3,
                                                                       19284, 5928, 19320, 398,
                                                                       408, 6492, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20472, 0, 3,
                                                                       19320, 5946, 19356, 408,
                                                                       418, 6522, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20532, 0, 3,
                                                                       19356, 5964, 19392, 418,
                                                                       428, 6552, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20592, 0, 3,
                                                                       19392, 5982, 19428, 428,
                                                                       438, 6582, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20652, 0, 3,
                                                                       19464, 6018, 19500, 458,
                                                                       468, 6612, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20712, 0, 3,
                                                                       19500, 6036, 19536, 468,
                                                                       478, 6642, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20772, 0, 3,
                                                                       19536, 6054, 19572, 478,
                                                                       488, 6672, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20832, 0, 3,
                                                                       19572, 6072, 19608, 488,
                                                                       498, 6702, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20892, 0, 3,
                                                                       19608, 6090, 19644, 498,
                                                                       508, 6732, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20952, 0, 3,
                                                                       19644, 6108, 19680, 508,
                                                                       518, 6762, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 21012, 0, 3,
                                                                       19680, 6126, 19716, 518,
                                                                       528, 6792, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 21072, 0, 3,
                                                                       19716, 6144, 19752, 528,
                                                                       538, 6822, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 21132, 0, 3,
                                                                       19752, 6162, 19788, 538,
                                                                       548, 6852, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 21192, 0, 3,
                                                                       19788, 6180, 19824, 548,
                                                                       558, 6882, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 21252, 0, 3,
                                                                       19824, 6198, 19860, 558,
                                                                       568, 6912, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 21312, 0, 3,
                                                                       19860, 6216, 19896, 568,
                                                                       578, 6942, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21372, 0, 3,
                                                                       19932, 6252, 19992, 598,
                                                                       613, 6972, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21462, 0, 3,
                                                                       19992, 6282, 20052, 613,
                                                                       628, 7017, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21552, 0, 3,
                                                                       20052, 6312, 20112, 628,
                                                                       643, 7062, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21642, 0, 3,
                                                                       20112, 6342, 20172, 643,
                                                                       658, 7107, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21732, 0, 3,
                                                                       20172, 6372, 20232, 658,
                                                                       673, 7152, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21822, 0, 3,
                                                                       20232, 6402, 20292, 673,
                                                                       688, 7197, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21912, 0, 3,
                                                                       20292, 6432, 20352, 688,
                                                                       703, 7242, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22002, 0, 3,
                                                                       20352, 6462, 20412, 703,
                                                                       718, 7287, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22092, 0, 3,
                                                                       20412, 6492, 20472, 718,
                                                                       733, 7332, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22182, 0, 3,
                                                                       20472, 6522, 20532, 733,
                                                                       748, 7377, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22272, 0, 3,
                                                                       20532, 6552, 20592, 748,
                                                                       763, 7422, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22362, 0, 3,
                                                                       20652, 6612, 20712, 793,
                                                                       808, 7467, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22452, 0, 3,
                                                                       20712, 6642, 20772, 808,
                                                                       823, 7512, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22542, 0, 3,
                                                                       20772, 6672, 20832, 823,
                                                                       838, 7557, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22632, 0, 3,
                                                                       20832, 6702, 20892, 838,
                                                                       853, 7602, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22722, 0, 3,
                                                                       20892, 6732, 20952, 853,
                                                                       868, 7647, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22812, 0, 3,
                                                                       20952, 6762, 21012, 868,
                                                                       883, 7692, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22902, 0, 3,
                                                                       21012, 6792, 21072, 883,
                                                                       898, 7737, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22992, 0, 3,
                                                                       21072, 6822, 21132, 898,
                                                                       913, 7782, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 23082, 0, 3,
                                                                       21132, 6852, 21192, 913,
                                                                       928, 7827, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 23172, 0, 3,
                                                                       21192, 6882, 21252, 928,
                                                                       943, 7872, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 23262, 0, 3,
                                                                       21252, 6912, 21312, 943,
                                                                       958, 7917, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23352, 0, 3,
                                                                       21372, 6972, 21462, 988,
                                                                       1009, 7962, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23478, 0, 3,
                                                                       21462, 7017, 21552, 1009,
                                                                       1030, 8025, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23604, 0, 3,
                                                                       21552, 7062, 21642, 1030,
                                                                       1051, 8088, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23730, 0, 3,
                                                                       21642, 7107, 21732, 1051,
                                                                       1072, 8151, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23856, 0, 3,
                                                                       21732, 7152, 21822, 1072,
                                                                       1093, 8214, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23982, 0, 3,
                                                                       21822, 7197, 21912, 1093,
                                                                       1114, 8277, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24108, 0, 3,
                                                                       21912, 7242, 22002, 1114,
                                                                       1135, 8340, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24234, 0, 3,
                                                                       22002, 7287, 22092, 1135,
                                                                       1156, 8403, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24360, 0, 3,
                                                                       22092, 7332, 22182, 1156,
                                                                       1177, 8466, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24486, 0, 3,
                                                                       22182, 7377, 22272, 1177,
                                                                       1198, 8529, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24612, 0, 3,
                                                                       22362, 7467, 22452, 1240,
                                                                       1261, 8592, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24738, 0, 3,
                                                                       22452, 7512, 22542, 1261,
                                                                       1282, 8655, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24864, 0, 3,
                                                                       22542, 7557, 22632, 1282,
                                                                       1303, 8718, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24990, 0, 3,
                                                                       22632, 7602, 22722, 1303,
                                                                       1324, 8781, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 25116, 0, 3,
                                                                       22722, 7647, 22812, 1324,
                                                                       1345, 8844, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 25242, 0, 3,
                                                                       22812, 7692, 22902, 1345,
                                                                       1366, 8907, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 25368, 0, 3,
                                                                       22902, 7737, 22992, 1366,
                                                                       1387, 8970, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 25494, 0, 3,
                                                                       22992, 7782, 23082, 1387,
                                                                       1408, 9033, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 25620, 0, 3,
                                                                       23082, 7827, 23172, 1408,
                                                                       1429, 9096, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 25746, 0, 3,
                                                                       23172, 7872, 23262, 1429,
                                                                       1450, 9159, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25872, 0, 3,
                                                                       23352, 7962, 23478, 1492,
                                                                       1520, 9222, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26040, 0, 3,
                                                                       23478, 8025, 23604, 1520,
                                                                       1548, 9306, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26208, 0, 3,
                                                                       23604, 8088, 23730, 1548,
                                                                       1576, 9390, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26376, 0, 3,
                                                                       23730, 8151, 23856, 1576,
                                                                       1604, 9474, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26544, 0, 3,
                                                                       23856, 8214, 23982, 1604,
                                                                       1632, 9558, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26712, 0, 3,
                                                                       23982, 8277, 24108, 1632,
                                                                       1660, 9642, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26880, 0, 3,
                                                                       24108, 8340, 24234, 1660,
                                                                       1688, 9726, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27048, 0, 3,
                                                                       24234, 8403, 24360, 1688,
                                                                       1716, 9810, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27216, 0, 3,
                                                                       24360, 8466, 24486, 1716,
                                                                       1744, 9894, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27384, 0, 3,
                                                                       24612, 8592, 24738, 1800,
                                                                       1828, 9978, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27552, 0, 3,
                                                                       24738, 8655, 24864, 1828,
                                                                       1856, 10062, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27720, 0, 3,
                                                                       24864, 8718, 24990, 1856,
                                                                       1884, 10146, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27888, 0, 3,
                                                                       24990, 8781, 25116, 1884,
                                                                       1912, 10230, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28056, 0, 3,
                                                                       25116, 8844, 25242, 1912,
                                                                       1940, 10314, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28224, 0, 3,
                                                                       25242, 8907, 25368, 1940,
                                                                       1968, 10398, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28392, 0, 3,
                                                                       25368, 8970, 25494, 1968,
                                                                       1996, 10482, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28560, 0, 3,
                                                                       25494, 9033, 25620, 1996,
                                                                       2024, 10566, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28728, 0, 3,
                                                                       25620, 9096, 25746, 2024,
                                                                       2052, 10650, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28896, 0, 3,
                                                                       25872, 9222, 26040, 2108,
                                                                       2144, 10734, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29112, 0, 3,
                                                                       26040, 9306, 26208, 2144,
                                                                       2180, 10842, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29328, 0, 3,
                                                                       26208, 9390, 26376, 2180,
                                                                       2216, 10950, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29544, 0, 3,
                                                                       26376, 9474, 26544, 2216,
                                                                       2252, 11058, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29760, 0, 3,
                                                                       26544, 9558, 26712, 2252,
                                                                       2288, 11166, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29976, 0, 3,
                                                                       26712, 9642, 26880, 2288,
                                                                       2324, 11274, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30192, 0, 3,
                                                                       26880, 9726, 27048, 2324,
                                                                       2360, 11382, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30408, 0, 3,
                                                                       27048, 9810, 27216, 2360,
                                                                       2396, 11490, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30624, 0, 3,
                                                                       27384, 9978, 27552, 2468,
                                                                       2504, 11598, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30840, 0, 3,
                                                                       27552, 10062, 27720, 2504,
                                                                       2540, 11706, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31056, 0, 3,
                                                                       27720, 10146, 27888, 2540,
                                                                       2576, 11814, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31272, 0, 3,
                                                                       27888, 10230, 28056, 2576,
                                                                       2612, 11922, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31488, 0, 3,
                                                                       28056, 10314, 28224, 2612,
                                                                       2648, 12030, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31704, 0, 3,
                                                                       28224, 10398, 28392, 2648,
                                                                       2684, 12138, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31920, 0, 3,
                                                                       28392, 10482, 28560, 2684,
                                                                       2720, 12246, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32136, 0, 3,
                                                                       28560, 10566, 28728, 2720,
                                                                       2756, 12354, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 32352, 0, 3,
                                                                       28896, 10734, 29112, 2828,
                                                                       2873, 12462, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 32622, 0, 3,
                                                                       29112, 10842, 29328, 2873,
                                                                       2918, 12597, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 32892, 0, 3,
                                                                       29328, 10950, 29544, 2918,
                                                                       2963, 12732, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33162, 0, 3,
                                                                       29544, 11058, 29760, 2963,
                                                                       3008, 12867, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33432, 0, 3,
                                                                       29760, 11166, 29976, 3008,
                                                                       3053, 13002, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33702, 0, 3,
                                                                       29976, 11274, 30192, 3053,
                                                                       3098, 13137, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33972, 0, 3,
                                                                       30192, 11382, 30408, 3098,
                                                                       3143, 13272, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 34242, 0, 3,
                                                                       30624, 11598, 30840, 3233,
                                                                       3278, 13407, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 34512, 0, 3,
                                                                       30840, 11706, 31056, 3278,
                                                                       3323, 13542, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 34782, 0, 3,
                                                                       31056, 11814, 31272, 3323,
                                                                       3368, 13677, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35052, 0, 3,
                                                                       31272, 11922, 31488, 3368,
                                                                       3413, 13812, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35322, 0, 3,
                                                                       31488, 12030, 31704, 3413,
                                                                       3458, 13947, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35592, 0, 3,
                                                                       31704, 12138, 31920, 3458,
                                                                       3503, 14082, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35862, 0, 3,
                                                                       31920, 12246, 32136, 3503,
                                                                       3548, 14217, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 36132, 0, 3,
                                                                       32352, 12462, 32622, 3638,
                                                                       3693, 14352, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 36462, 0, 3,
                                                                       32622, 12597, 32892, 3693,
                                                                       3748, 14517, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 36792, 0, 3,
                                                                       32892, 12732, 33162, 3748,
                                                                       3803, 14682, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 37122, 0, 3,
                                                                       33162, 12867, 33432, 3803,
                                                                       3858, 14847, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 37452, 0, 3,
                                                                       33432, 13002, 33702, 3858,
                                                                       3913, 15012, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 37782, 0, 3,
                                                                       33702, 13137, 33972, 3913,
                                                                       3968, 15177, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 38112, 0, 3,
                                                                       34242, 13407, 34512, 4078,
                                                                       4133, 15342, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 38442, 0, 3,
                                                                       34512, 13542, 34782, 4133,
                                                                       4188, 15507, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 38772, 0, 3,
                                                                       34782, 13677, 35052, 4188,
                                                                       4243, 15672, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 39102, 0, 3,
                                                                       35052, 13812, 35322, 4243,
                                                                       4298, 15837, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 39432, 0, 3,
                                                                       35322, 13947, 35592, 4298,
                                                                       4353, 16002, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 39762, 0, 3,
                                                                       35592, 14082, 35862, 4353,
                                                                       4408, 16167, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 40092, 0, 3,
                                                                       36132, 14352, 36462, 4518,
                                                                       4584, 16332, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 40488, 0, 3,
                                                                       36462, 14517, 36792, 4584,
                                                                       4650, 16530, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 40884, 0, 3,
                                                                       36792, 14682, 37122, 4650,
                                                                       4716, 16728, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 41280, 0, 3,
                                                                       37122, 14847, 37452, 4716,
                                                                       4782, 16926, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 41676, 0, 3,
                                                                       37452, 15012, 37782, 4782,
                                                                       4848, 17124, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 42072, 0, 3,
                                                                       38112, 15342, 38442, 4980,
                                                                       5046, 17322, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 42468, 0, 3,
                                                                       38442, 15507, 38772, 5046,
                                                                       5112, 17520, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 42864, 0, 3,
                                                                       38772, 15672, 39102, 5112,
                                                                       5178, 17718, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 43260, 0, 3,
                                                                       39102, 15837, 39432, 5178,
                                                                       5244, 17916, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 43656, 0, 3,
                                                                       39432, 16002, 39762, 5244,
                                                                       5310, 18114, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44052, 3, 5442,
                                                                       5445, 18324, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44062, 3, 5445,
                                                                       5448, 18330, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44072, 3, 5448,
                                                                       5451, 18336, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44082, 3, 5451,
                                                                       5454, 18342, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44092, 3, 5454,
                                                                       5457, 18348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44102, 3, 5457,
                                                                       5460, 18354, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44112, 3, 5460,
                                                                       5463, 18360, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44122, 3, 5463,
                                                                       5466, 18366, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44132, 3, 5466,
                                                                       5469, 18372, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44142, 3, 5469,
                                                                       5472, 18378, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44152, 3, 5472,
                                                                       5475, 18384, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44162, 3, 5475,
                                                                       5478, 18390, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44172, 3, 5478,
                                                                       5481, 18396, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44182, 3, 5487,
                                                                       5490, 18414, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44192, 3, 5490,
                                                                       5493, 18420, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44202, 3, 5493,
                                                                       5496, 18426, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44212, 3, 5496,
                                                                       5499, 18432, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44222, 3, 5499,
                                                                       5502, 18438, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44232, 3, 5502,
                                                                       5505, 18444, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44242, 3, 5505,
                                                                       5508, 18450, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44252, 3, 5508,
                                                                       5511, 18456, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44262, 3, 5511,
                                                                       5514, 18462, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44272, 3, 5514,
                                                                       5517, 18468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44282, 3, 5517,
                                                                       5520, 18474, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44292, 3, 5520,
                                                                       5523, 18480, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 44302, 3, 5523,
                                                                       5526, 18486, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44312, 0, 3,
                                                                       44052, 18324, 44062, 5532,
                                                                       5541, 18528, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44342, 0, 3,
                                                                       44062, 18330, 44072, 5541,
                                                                       5550, 18546, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44372, 0, 3,
                                                                       44072, 18336, 44082, 5550,
                                                                       5559, 18564, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44402, 0, 3,
                                                                       44082, 18342, 44092, 5559,
                                                                       5568, 18582, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44432, 0, 3,
                                                                       44092, 18348, 44102, 5568,
                                                                       5577, 18600, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44462, 0, 3,
                                                                       44102, 18354, 44112, 5577,
                                                                       5586, 18618, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44492, 0, 3,
                                                                       44112, 18360, 44122, 5586,
                                                                       5595, 18636, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44522, 0, 3,
                                                                       44122, 18366, 44132, 5595,
                                                                       5604, 18654, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44552, 0, 3,
                                                                       44132, 18372, 44142, 5604,
                                                                       5613, 18672, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44582, 0, 3,
                                                                       44142, 18378, 44152, 5613,
                                                                       5622, 18690, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44612, 0, 3,
                                                                       44152, 18384, 44162, 5622,
                                                                       5631, 18708, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44642, 0, 3,
                                                                       44162, 18390, 44172, 5631,
                                                                       5640, 18726, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44672, 0, 3,
                                                                       44182, 18414, 44192, 5658,
                                                                       5667, 18780, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44702, 0, 3,
                                                                       44192, 18420, 44202, 5667,
                                                                       5676, 18798, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44732, 0, 3,
                                                                       44202, 18426, 44212, 5676,
                                                                       5685, 18816, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44762, 0, 3,
                                                                       44212, 18432, 44222, 5685,
                                                                       5694, 18834, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44792, 0, 3,
                                                                       44222, 18438, 44232, 5694,
                                                                       5703, 18852, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44822, 0, 3,
                                                                       44232, 18444, 44242, 5703,
                                                                       5712, 18870, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44852, 0, 3,
                                                                       44242, 18450, 44252, 5712,
                                                                       5721, 18888, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44882, 0, 3,
                                                                       44252, 18456, 44262, 5721,
                                                                       5730, 18906, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44912, 0, 3,
                                                                       44262, 18462, 44272, 5730,
                                                                       5739, 18924, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44942, 0, 3,
                                                                       44272, 18468, 44282, 5739,
                                                                       5748, 18942, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44972, 0, 3,
                                                                       44282, 18474, 44292, 5748,
                                                                       5757, 18960, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 45002, 0, 3,
                                                                       44292, 18480, 44302, 5757,
                                                                       5766, 18978, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45032, 0, 3,
                                                                       44312, 18528, 44342, 5784,
                                                                       5802, 19068, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45092, 0, 3,
                                                                       44342, 18546, 44372, 5802,
                                                                       5820, 19104, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45152, 0, 3,
                                                                       44372, 18564, 44402, 5820,
                                                                       5838, 19140, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45212, 0, 3,
                                                                       44402, 18582, 44432, 5838,
                                                                       5856, 19176, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45272, 0, 3,
                                                                       44432, 18600, 44462, 5856,
                                                                       5874, 19212, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45332, 0, 3,
                                                                       44462, 18618, 44492, 5874,
                                                                       5892, 19248, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45392, 0, 3,
                                                                       44492, 18636, 44522, 5892,
                                                                       5910, 19284, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45452, 0, 3,
                                                                       44522, 18654, 44552, 5910,
                                                                       5928, 19320, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45512, 0, 3,
                                                                       44552, 18672, 44582, 5928,
                                                                       5946, 19356, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45572, 0, 3,
                                                                       44582, 18690, 44612, 5946,
                                                                       5964, 19392, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45632, 0, 3,
                                                                       44612, 18708, 44642, 5964,
                                                                       5982, 19428, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45692, 0, 3,
                                                                       44672, 18780, 44702, 6018,
                                                                       6036, 19536, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45752, 0, 3,
                                                                       44702, 18798, 44732, 6036,
                                                                       6054, 19572, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45812, 0, 3,
                                                                       44732, 18816, 44762, 6054,
                                                                       6072, 19608, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45872, 0, 3,
                                                                       44762, 18834, 44792, 6072,
                                                                       6090, 19644, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45932, 0, 3,
                                                                       44792, 18852, 44822, 6090,
                                                                       6108, 19680, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45992, 0, 3,
                                                                       44822, 18870, 44852, 6108,
                                                                       6126, 19716, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 46052, 0, 3,
                                                                       44852, 18888, 44882, 6126,
                                                                       6144, 19752, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 46112, 0, 3,
                                                                       44882, 18906, 44912, 6144,
                                                                       6162, 19788, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 46172, 0, 3,
                                                                       44912, 18924, 44942, 6162,
                                                                       6180, 19824, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 46232, 0, 3,
                                                                       44942, 18942, 44972, 6180,
                                                                       6198, 19860, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 46292, 0, 3,
                                                                       44972, 18960, 45002, 6198,
                                                                       6216, 19896, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46352, 0, 3,
                                                                       45032, 19068, 45092, 6252,
                                                                       6282, 20052, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46452, 0, 3,
                                                                       45092, 19104, 45152, 6282,
                                                                       6312, 20112, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46552, 0, 3,
                                                                       45152, 19140, 45212, 6312,
                                                                       6342, 20172, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46652, 0, 3,
                                                                       45212, 19176, 45272, 6342,
                                                                       6372, 20232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46752, 0, 3,
                                                                       45272, 19212, 45332, 6372,
                                                                       6402, 20292, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46852, 0, 3,
                                                                       45332, 19248, 45392, 6402,
                                                                       6432, 20352, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46952, 0, 3,
                                                                       45392, 19284, 45452, 6432,
                                                                       6462, 20412, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47052, 0, 3,
                                                                       45452, 19320, 45512, 6462,
                                                                       6492, 20472, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47152, 0, 3,
                                                                       45512, 19356, 45572, 6492,
                                                                       6522, 20532, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47252, 0, 3,
                                                                       45572, 19392, 45632, 6522,
                                                                       6552, 20592, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47352, 0, 3,
                                                                       45692, 19536, 45752, 6612,
                                                                       6642, 20772, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47452, 0, 3,
                                                                       45752, 19572, 45812, 6642,
                                                                       6672, 20832, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47552, 0, 3,
                                                                       45812, 19608, 45872, 6672,
                                                                       6702, 20892, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47652, 0, 3,
                                                                       45872, 19644, 45932, 6702,
                                                                       6732, 20952, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47752, 0, 3,
                                                                       45932, 19680, 45992, 6732,
                                                                       6762, 21012, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47852, 0, 3,
                                                                       45992, 19716, 46052, 6762,
                                                                       6792, 21072, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47952, 0, 3,
                                                                       46052, 19752, 46112, 6792,
                                                                       6822, 21132, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48052, 0, 3,
                                                                       46112, 19788, 46172, 6822,
                                                                       6852, 21192, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48152, 0, 3,
                                                                       46172, 19824, 46232, 6852,
                                                                       6882, 21252, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48252, 0, 3,
                                                                       46232, 19860, 46292, 6882,
                                                                       6912, 21312, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48352, 0, 3,
                                                                       46352, 20052, 46452, 6972,
                                                                       7017, 21552, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48502, 0, 3,
                                                                       46452, 20112, 46552, 7017,
                                                                       7062, 21642, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48652, 0, 3,
                                                                       46552, 20172, 46652, 7062,
                                                                       7107, 21732, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48802, 0, 3,
                                                                       46652, 20232, 46752, 7107,
                                                                       7152, 21822, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48952, 0, 3,
                                                                       46752, 20292, 46852, 7152,
                                                                       7197, 21912, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49102, 0, 3,
                                                                       46852, 20352, 46952, 7197,
                                                                       7242, 22002, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49252, 0, 3,
                                                                       46952, 20412, 47052, 7242,
                                                                       7287, 22092, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49402, 0, 3,
                                                                       47052, 20472, 47152, 7287,
                                                                       7332, 22182, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49552, 0, 3,
                                                                       47152, 20532, 47252, 7332,
                                                                       7377, 22272, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49702, 0, 3,
                                                                       47352, 20772, 47452, 7467,
                                                                       7512, 22542, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49852, 0, 3,
                                                                       47452, 20832, 47552, 7512,
                                                                       7557, 22632, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50002, 0, 3,
                                                                       47552, 20892, 47652, 7557,
                                                                       7602, 22722, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50152, 0, 3,
                                                                       47652, 20952, 47752, 7602,
                                                                       7647, 22812, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50302, 0, 3,
                                                                       47752, 21012, 47852, 7647,
                                                                       7692, 22902, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50452, 0, 3,
                                                                       47852, 21072, 47952, 7692,
                                                                       7737, 22992, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50602, 0, 3,
                                                                       47952, 21132, 48052, 7737,
                                                                       7782, 23082, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50752, 0, 3,
                                                                       48052, 21192, 48152, 7782,
                                                                       7827, 23172, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50902, 0, 3,
                                                                       48152, 21252, 48252, 7827,
                                                                       7872, 23262, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51052, 0, 3,
                                                                       48352, 21552, 48502, 7962,
                                                                       8025, 23604, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51262, 0, 3,
                                                                       48502, 21642, 48652, 8025,
                                                                       8088, 23730, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51472, 0, 3,
                                                                       48652, 21732, 48802, 8088,
                                                                       8151, 23856, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51682, 0, 3,
                                                                       48802, 21822, 48952, 8151,
                                                                       8214, 23982, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51892, 0, 3,
                                                                       48952, 21912, 49102, 8214,
                                                                       8277, 24108, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52102, 0, 3,
                                                                       49102, 22002, 49252, 8277,
                                                                       8340, 24234, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52312, 0, 3,
                                                                       49252, 22092, 49402, 8340,
                                                                       8403, 24360, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52522, 0, 3,
                                                                       49402, 22182, 49552, 8403,
                                                                       8466, 24486, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52732, 0, 3,
                                                                       49702, 22542, 49852, 8592,
                                                                       8655, 24864, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52942, 0, 3,
                                                                       49852, 22632, 50002, 8655,
                                                                       8718, 24990, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53152, 0, 3,
                                                                       50002, 22722, 50152, 8718,
                                                                       8781, 25116, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53362, 0, 3,
                                                                       50152, 22812, 50302, 8781,
                                                                       8844, 25242, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53572, 0, 3,
                                                                       50302, 22902, 50452, 8844,
                                                                       8907, 25368, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53782, 0, 3,
                                                                       50452, 22992, 50602, 8907,
                                                                       8970, 25494, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53992, 0, 3,
                                                                       50602, 23082, 50752, 8970,
                                                                       9033, 25620, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 54202, 0, 3,
                                                                       50752, 23172, 50902, 9033,
                                                                       9096, 25746, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54412, 0, 3,
                                                                       51052, 23604, 51262, 9222,
                                                                       9306, 26208, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54692, 0, 3,
                                                                       51262, 23730, 51472, 9306,
                                                                       9390, 26376, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54972, 0, 3,
                                                                       51472, 23856, 51682, 9390,
                                                                       9474, 26544, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55252, 0, 3,
                                                                       51682, 23982, 51892, 9474,
                                                                       9558, 26712, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55532, 0, 3,
                                                                       51892, 24108, 52102, 9558,
                                                                       9642, 26880, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55812, 0, 3,
                                                                       52102, 24234, 52312, 9642,
                                                                       9726, 27048, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56092, 0, 3,
                                                                       52312, 24360, 52522, 9726,
                                                                       9810, 27216, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56372, 0, 3,
                                                                       52732, 24864, 52942, 9978,
                                                                       10062, 27720, ncols,
                                                                       gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56652, 0, 3,
                                                                       52942, 24990, 53152,
                                                                       10062, 10146, 27888,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56932, 0, 3,
                                                                       53152, 25116, 53362,
                                                                       10146, 10230, 28056,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 57212, 0, 3,
                                                                       53362, 25242, 53572,
                                                                       10230, 10314, 28224,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 57492, 0, 3,
                                                                       53572, 25368, 53782,
                                                                       10314, 10398, 28392,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 57772, 0, 3,
                                                                       53782, 25494, 53992,
                                                                       10398, 10482, 28560,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 58052, 0, 3,
                                                                       53992, 25620, 54202,
                                                                       10482, 10566, 28728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58332, 0, 3,
                                                                       54412, 26208, 54692,
                                                                       10734, 10842, 29328,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58692, 0, 3,
                                                                       54692, 26376, 54972,
                                                                       10842, 10950, 29544,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59052, 0, 3,
                                                                       54972, 26544, 55252,
                                                                       10950, 11058, 29760,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59412, 0, 3,
                                                                       55252, 26712, 55532,
                                                                       11058, 11166, 29976,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59772, 0, 3,
                                                                       55532, 26880, 55812,
                                                                       11166, 11274, 30192,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 60132, 0, 3,
                                                                       55812, 27048, 56092,
                                                                       11274, 11382, 30408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 60492, 0, 3,
                                                                       56372, 27720, 56652,
                                                                       11598, 11706, 31056,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 60852, 0, 3,
                                                                       56652, 27888, 56932,
                                                                       11706, 11814, 31272,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 61212, 0, 3,
                                                                       56932, 28056, 57212,
                                                                       11814, 11922, 31488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 61572, 0, 3,
                                                                       57212, 28224, 57492,
                                                                       11922, 12030, 31704,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 61932, 0, 3,
                                                                       57492, 28392, 57772,
                                                                       12030, 12138, 31920,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 62292, 0, 3,
                                                                       57772, 28560, 58052,
                                                                       12138, 12246, 32136,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 62652, 0, 3,
                                                                       58332, 29328, 58692,
                                                                       12462, 12597, 32892,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 63102, 0, 3,
                                                                       58692, 29544, 59052,
                                                                       12597, 12732, 33162,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 63552, 0, 3,
                                                                       59052, 29760, 59412,
                                                                       12732, 12867, 33432,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 64002, 0, 3,
                                                                       59412, 29976, 59772,
                                                                       12867, 13002, 33702,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 64452, 0, 3,
                                                                       59772, 30192, 60132,
                                                                       13002, 13137, 33972,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 64902, 0, 3,
                                                                       60492, 31056, 60852,
                                                                       13407, 13542, 34782,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 65352, 0, 3,
                                                                       60852, 31272, 61212,
                                                                       13542, 13677, 35052,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 65802, 0, 3,
                                                                       61212, 31488, 61572,
                                                                       13677, 13812, 35322,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 66252, 0, 3,
                                                                       61572, 31704, 61932,
                                                                       13812, 13947, 35592,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 66702, 0, 3,
                                                                       61932, 31920, 62292,
                                                                       13947, 14082, 35862,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 67152, 0, 3,
                                                                       62652, 32892, 63102,
                                                                       14352, 14517, 36792,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 67702, 0, 3,
                                                                       63102, 33162, 63552,
                                                                       14517, 14682, 37122,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 68252, 0, 3,
                                                                       63552, 33432, 64002,
                                                                       14682, 14847, 37452,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 68802, 0, 3,
                                                                       64002, 33702, 64452,
                                                                       14847, 15012, 37782,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 69352, 0, 3,
                                                                       64902, 34782, 65352,
                                                                       15342, 15507, 38772,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 69902, 0, 3,
                                                                       65352, 35052, 65802,
                                                                       15507, 15672, 39102,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 70452, 0, 3,
                                                                       65802, 35322, 66252,
                                                                       15672, 15837, 39432,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 71002, 0, 3,
                                                                       66252, 35592, 66702,
                                                                       15837, 16002, 39762,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 71552, 0, 3,
                                                                       67152, 36792, 67702,
                                                                       16332, 16530, 40884,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 72212, 0, 3,
                                                                       67702, 37122, 68252,
                                                                       16530, 16728, 41280,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 72872, 0, 3,
                                                                       68252, 37452, 68802,
                                                                       16728, 16926, 41676,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 73532, 0, 3,
                                                                       69352, 38772, 69902,
                                                                       17322, 17520, 42864,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 74192, 0, 3,
                                                                       69902, 39102, 70452,
                                                                       17520, 17718, 43260,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 74852, 0, 3,
                                                                       70452, 39432, 71002,
                                                                       17718, 17916, 43656,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75512, 3, 18312,
                                                                       18318, 44052, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75527, 3, 18318,
                                                                       18324, 44062, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75542, 3, 18324,
                                                                       18330, 44072, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75557, 3, 18330,
                                                                       18336, 44082, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75572, 3, 18336,
                                                                       18342, 44092, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75587, 3, 18342,
                                                                       18348, 44102, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75602, 3, 18348,
                                                                       18354, 44112, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75617, 3, 18354,
                                                                       18360, 44122, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75632, 3, 18360,
                                                                       18366, 44132, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75647, 3, 18366,
                                                                       18372, 44142, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75662, 3, 18372,
                                                                       18378, 44152, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75677, 3, 18378,
                                                                       18384, 44162, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75692, 3, 18384,
                                                                       18390, 44172, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75707, 3, 18402,
                                                                       18408, 44182, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75722, 3, 18408,
                                                                       18414, 44192, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75737, 3, 18414,
                                                                       18420, 44202, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75752, 3, 18420,
                                                                       18426, 44212, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75767, 3, 18426,
                                                                       18432, 44222, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75782, 3, 18432,
                                                                       18438, 44232, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75797, 3, 18438,
                                                                       18444, 44242, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75812, 3, 18444,
                                                                       18450, 44252, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75827, 3, 18450,
                                                                       18456, 44262, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75842, 3, 18456,
                                                                       18462, 44272, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75857, 3, 18462,
                                                                       18468, 44282, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75872, 3, 18468,
                                                                       18474, 44292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 75887, 3, 18474,
                                                                       18480, 44302, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 75902, 0, 3,
                                                                       75512, 44052, 75527,
                                                                       18492, 18510, 44312,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 75947, 0, 3,
                                                                       75527, 44062, 75542,
                                                                       18510, 18528, 44342,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 75992, 0, 3,
                                                                       75542, 44072, 75557,
                                                                       18528, 18546, 44372,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76037, 0, 3,
                                                                       75557, 44082, 75572,
                                                                       18546, 18564, 44402,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76082, 0, 3,
                                                                       75572, 44092, 75587,
                                                                       18564, 18582, 44432,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76127, 0, 3,
                                                                       75587, 44102, 75602,
                                                                       18582, 18600, 44462,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76172, 0, 3,
                                                                       75602, 44112, 75617,
                                                                       18600, 18618, 44492,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76217, 0, 3,
                                                                       75617, 44122, 75632,
                                                                       18618, 18636, 44522,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76262, 0, 3,
                                                                       75632, 44132, 75647,
                                                                       18636, 18654, 44552,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76307, 0, 3,
                                                                       75647, 44142, 75662,
                                                                       18654, 18672, 44582,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76352, 0, 3,
                                                                       75662, 44152, 75677,
                                                                       18672, 18690, 44612,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76397, 0, 3,
                                                                       75677, 44162, 75692,
                                                                       18690, 18708, 44642,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76442, 0, 3,
                                                                       75707, 44182, 75722,
                                                                       18744, 18762, 44672,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76487, 0, 3,
                                                                       75722, 44192, 75737,
                                                                       18762, 18780, 44702,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76532, 0, 3,
                                                                       75737, 44202, 75752,
                                                                       18780, 18798, 44732,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76577, 0, 3,
                                                                       75752, 44212, 75767,
                                                                       18798, 18816, 44762,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76622, 0, 3,
                                                                       75767, 44222, 75782,
                                                                       18816, 18834, 44792,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76667, 0, 3,
                                                                       75782, 44232, 75797,
                                                                       18834, 18852, 44822,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76712, 0, 3,
                                                                       75797, 44242, 75812,
                                                                       18852, 18870, 44852,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76757, 0, 3,
                                                                       75812, 44252, 75827,
                                                                       18870, 18888, 44882,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76802, 0, 3,
                                                                       75827, 44262, 75842,
                                                                       18888, 18906, 44912,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76847, 0, 3,
                                                                       75842, 44272, 75857,
                                                                       18906, 18924, 44942,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76892, 0, 3,
                                                                       75857, 44282, 75872,
                                                                       18924, 18942, 44972,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76937, 0, 3,
                                                                       75872, 44292, 75887,
                                                                       18942, 18960, 45002,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76982, 0, 3,
                                                                       75902, 44312, 75947,
                                                                       18996, 19032, 45032,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77072, 0, 3,
                                                                       75947, 44342, 75992,
                                                                       19032, 19068, 45092,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77162, 0, 3,
                                                                       75992, 44372, 76037,
                                                                       19068, 19104, 45152,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77252, 0, 3,
                                                                       76037, 44402, 76082,
                                                                       19104, 19140, 45212,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77342, 0, 3,
                                                                       76082, 44432, 76127,
                                                                       19140, 19176, 45272,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77432, 0, 3,
                                                                       76127, 44462, 76172,
                                                                       19176, 19212, 45332,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77522, 0, 3,
                                                                       76172, 44492, 76217,
                                                                       19212, 19248, 45392,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77612, 0, 3,
                                                                       76217, 44522, 76262,
                                                                       19248, 19284, 45452,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77702, 0, 3,
                                                                       76262, 44552, 76307,
                                                                       19284, 19320, 45512,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77792, 0, 3,
                                                                       76307, 44582, 76352,
                                                                       19320, 19356, 45572,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77882, 0, 3,
                                                                       76352, 44612, 76397,
                                                                       19356, 19392, 45632,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77972, 0, 3,
                                                                       76442, 44672, 76487,
                                                                       19464, 19500, 45692,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78062, 0, 3,
                                                                       76487, 44702, 76532,
                                                                       19500, 19536, 45752,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78152, 0, 3,
                                                                       76532, 44732, 76577,
                                                                       19536, 19572, 45812,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78242, 0, 3,
                                                                       76577, 44762, 76622,
                                                                       19572, 19608, 45872,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78332, 0, 3,
                                                                       76622, 44792, 76667,
                                                                       19608, 19644, 45932,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78422, 0, 3,
                                                                       76667, 44822, 76712,
                                                                       19644, 19680, 45992,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78512, 0, 3,
                                                                       76712, 44852, 76757,
                                                                       19680, 19716, 46052,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78602, 0, 3,
                                                                       76757, 44882, 76802,
                                                                       19716, 19752, 46112,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78692, 0, 3,
                                                                       76802, 44912, 76847,
                                                                       19752, 19788, 46172,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78782, 0, 3,
                                                                       76847, 44942, 76892,
                                                                       19788, 19824, 46232,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78872, 0, 3,
                                                                       76892, 44972, 76937,
                                                                       19824, 19860, 46292,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78962, 0, 3,
                                                                       76982, 45032, 77072,
                                                                       19932, 19992, 46352,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79112, 0, 3,
                                                                       77072, 45092, 77162,
                                                                       19992, 20052, 46452,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79262, 0, 3,
                                                                       77162, 45152, 77252,
                                                                       20052, 20112, 46552,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79412, 0, 3,
                                                                       77252, 45212, 77342,
                                                                       20112, 20172, 46652,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79562, 0, 3,
                                                                       77342, 45272, 77432,
                                                                       20172, 20232, 46752,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79712, 0, 3,
                                                                       77432, 45332, 77522,
                                                                       20232, 20292, 46852,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79862, 0, 3,
                                                                       77522, 45392, 77612,
                                                                       20292, 20352, 46952,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80012, 0, 3,
                                                                       77612, 45452, 77702,
                                                                       20352, 20412, 47052,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80162, 0, 3,
                                                                       77702, 45512, 77792,
                                                                       20412, 20472, 47152,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80312, 0, 3,
                                                                       77792, 45572, 77882,
                                                                       20472, 20532, 47252,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80462, 0, 3,
                                                                       77972, 45692, 78062,
                                                                       20652, 20712, 47352,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80612, 0, 3,
                                                                       78062, 45752, 78152,
                                                                       20712, 20772, 47452,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80762, 0, 3,
                                                                       78152, 45812, 78242,
                                                                       20772, 20832, 47552,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80912, 0, 3,
                                                                       78242, 45872, 78332,
                                                                       20832, 20892, 47652,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 81062, 0, 3,
                                                                       78332, 45932, 78422,
                                                                       20892, 20952, 47752,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 81212, 0, 3,
                                                                       78422, 45992, 78512,
                                                                       20952, 21012, 47852,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 81362, 0, 3,
                                                                       78512, 46052, 78602,
                                                                       21012, 21072, 47952,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 81512, 0, 3,
                                                                       78602, 46112, 78692,
                                                                       21072, 21132, 48052,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 81662, 0, 3,
                                                                       78692, 46172, 78782,
                                                                       21132, 21192, 48152,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 81812, 0, 3,
                                                                       78782, 46232, 78872,
                                                                       21192, 21252, 48252,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81962, 0, 3,
                                                                       78962, 46352, 79112,
                                                                       21372, 21462, 48352,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82187, 0, 3,
                                                                       79112, 46452, 79262,
                                                                       21462, 21552, 48502,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82412, 0, 3,
                                                                       79262, 46552, 79412,
                                                                       21552, 21642, 48652,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82637, 0, 3,
                                                                       79412, 46652, 79562,
                                                                       21642, 21732, 48802,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82862, 0, 3,
                                                                       79562, 46752, 79712,
                                                                       21732, 21822, 48952,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 83087, 0, 3,
                                                                       79712, 46852, 79862,
                                                                       21822, 21912, 49102,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 83312, 0, 3,
                                                                       79862, 46952, 80012,
                                                                       21912, 22002, 49252,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 83537, 0, 3,
                                                                       80012, 47052, 80162,
                                                                       22002, 22092, 49402,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 83762, 0, 3,
                                                                       80162, 47152, 80312,
                                                                       22092, 22182, 49552,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 83987, 0, 3,
                                                                       80462, 47352, 80612,
                                                                       22362, 22452, 49702,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 84212, 0, 3,
                                                                       80612, 47452, 80762,
                                                                       22452, 22542, 49852,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 84437, 0, 3,
                                                                       80762, 47552, 80912,
                                                                       22542, 22632, 50002,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 84662, 0, 3,
                                                                       80912, 47652, 81062,
                                                                       22632, 22722, 50152,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 84887, 0, 3,
                                                                       81062, 47752, 81212,
                                                                       22722, 22812, 50302,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 85112, 0, 3,
                                                                       81212, 47852, 81362,
                                                                       22812, 22902, 50452,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 85337, 0, 3,
                                                                       81362, 47952, 81512,
                                                                       22902, 22992, 50602,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 85562, 0, 3,
                                                                       81512, 48052, 81662,
                                                                       22992, 23082, 50752,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 85787, 0, 3,
                                                                       81662, 48152, 81812,
                                                                       23082, 23172, 50902,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 86012, 0, 3,
                                                                       81962, 48352, 82187,
                                                                       23352, 23478, 51052,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 86327, 0, 3,
                                                                       82187, 48502, 82412,
                                                                       23478, 23604, 51262,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 86642, 0, 3,
                                                                       82412, 48652, 82637,
                                                                       23604, 23730, 51472,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 86957, 0, 3,
                                                                       82637, 48802, 82862,
                                                                       23730, 23856, 51682,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 87272, 0, 3,
                                                                       82862, 48952, 83087,
                                                                       23856, 23982, 51892,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 87587, 0, 3,
                                                                       83087, 49102, 83312,
                                                                       23982, 24108, 52102,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 87902, 0, 3,
                                                                       83312, 49252, 83537,
                                                                       24108, 24234, 52312,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 88217, 0, 3,
                                                                       83537, 49402, 83762,
                                                                       24234, 24360, 52522,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 88532, 0, 3,
                                                                       83987, 49702, 84212,
                                                                       24612, 24738, 52732,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 88847, 0, 3,
                                                                       84212, 49852, 84437,
                                                                       24738, 24864, 52942,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 89162, 0, 3,
                                                                       84437, 50002, 84662,
                                                                       24864, 24990, 53152,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 89477, 0, 3,
                                                                       84662, 50152, 84887,
                                                                       24990, 25116, 53362,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 89792, 0, 3,
                                                                       84887, 50302, 85112,
                                                                       25116, 25242, 53572,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 90107, 0, 3,
                                                                       85112, 50452, 85337,
                                                                       25242, 25368, 53782,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 90422, 0, 3,
                                                                       85337, 50602, 85562,
                                                                       25368, 25494, 53992,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 90737, 0, 3,
                                                                       85562, 50752, 85787,
                                                                       25494, 25620, 54202,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 91052, 0, 3,
                                                                       86012, 51052, 86327,
                                                                       25872, 26040, 54412,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 91472, 0, 3,
                                                                       86327, 51262, 86642,
                                                                       26040, 26208, 54692,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 91892, 0, 3,
                                                                       86642, 51472, 86957,
                                                                       26208, 26376, 54972,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 92312, 0, 3,
                                                                       86957, 51682, 87272,
                                                                       26376, 26544, 55252,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 92732, 0, 3,
                                                                       87272, 51892, 87587,
                                                                       26544, 26712, 55532,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 93152, 0, 3,
                                                                       87587, 52102, 87902,
                                                                       26712, 26880, 55812,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 93572, 0, 3,
                                                                       87902, 52312, 88217,
                                                                       26880, 27048, 56092,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 93992, 0, 3,
                                                                       88532, 52732, 88847,
                                                                       27384, 27552, 56372,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 94412, 0, 3,
                                                                       88847, 52942, 89162,
                                                                       27552, 27720, 56652,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 94832, 0, 3,
                                                                       89162, 53152, 89477,
                                                                       27720, 27888, 56932,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 95252, 0, 3,
                                                                       89477, 53362, 89792,
                                                                       27888, 28056, 57212,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 95672, 0, 3,
                                                                       89792, 53572, 90107,
                                                                       28056, 28224, 57492,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 96092, 0, 3,
                                                                       90107, 53782, 90422,
                                                                       28224, 28392, 57772,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 96512, 0, 3,
                                                                       90422, 53992, 90737,
                                                                       28392, 28560, 58052,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 96932, 0, 3,
                                                                       91052, 54412, 91472,
                                                                       28896, 29112, 58332,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 97472, 0, 3,
                                                                       91472, 54692, 91892,
                                                                       29112, 29328, 58692,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 98012, 0, 3,
                                                                       91892, 54972, 92312,
                                                                       29328, 29544, 59052,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 98552, 0, 3,
                                                                       92312, 55252, 92732,
                                                                       29544, 29760, 59412,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 99092, 0, 3,
                                                                       92732, 55532, 93152,
                                                                       29760, 29976, 59772,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 99632, 0, 3,
                                                                       93152, 55812, 93572,
                                                                       29976, 30192, 60132,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 100172, 0, 3,
                                                                       93992, 56372, 94412,
                                                                       30624, 30840, 60492,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 100712, 0, 3,
                                                                       94412, 56652, 94832,
                                                                       30840, 31056, 60852,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 101252, 0, 3,
                                                                       94832, 56932, 95252,
                                                                       31056, 31272, 61212,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 101792, 0, 3,
                                                                       95252, 57212, 95672,
                                                                       31272, 31488, 61572,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 102332, 0, 3,
                                                                       95672, 57492, 96092,
                                                                       31488, 31704, 61932,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 102872, 0, 3,
                                                                       96092, 57772, 96512,
                                                                       31704, 31920, 62292,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 103412, 0, 3,
                                                                       96932, 58332, 97472,
                                                                       32352, 32622, 62652,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 104087, 0, 3,
                                                                       97472, 58692, 98012,
                                                                       32622, 32892, 63102,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 104762, 0, 3,
                                                                       98012, 59052, 98552,
                                                                       32892, 33162, 63552,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 105437, 0, 3,
                                                                       98552, 59412, 99092,
                                                                       33162, 33432, 64002,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 106112, 0, 3,
                                                                       99092, 59772, 99632,
                                                                       33432, 33702, 64452,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 106787, 0, 3,
                                                                       100172, 60492, 100712,
                                                                       34242, 34512, 64902,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 107462, 0, 3,
                                                                       100712, 60852, 101252,
                                                                       34512, 34782, 65352,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 108137, 0, 3,
                                                                       101252, 61212, 101792,
                                                                       34782, 35052, 65802,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 108812, 0, 3,
                                                                       101792, 61572, 102332,
                                                                       35052, 35322, 66252,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 109487, 0, 3,
                                                                       102332, 61932, 102872,
                                                                       35322, 35592, 66702,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 110162, 0, 3,
                                                                       103412, 62652, 104087,
                                                                       36132, 36462, 67152,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 110987, 0, 3,
                                                                       104087, 63102, 104762,
                                                                       36462, 36792, 67702,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 111812, 0, 3,
                                                                       104762, 63552, 105437,
                                                                       36792, 37122, 68252,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 112637, 0, 3,
                                                                       105437, 64002, 106112,
                                                                       37122, 37452, 68802,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 113462, 0, 3,
                                                                       106787, 64902, 107462,
                                                                       38112, 38442, 69352,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 114287, 0, 3,
                                                                       107462, 65352, 108137,
                                                                       38442, 38772, 69902,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 115112, 0, 3,
                                                                       108137, 65802, 108812,
                                                                       38772, 39102, 70452,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 115937, 0, 3,
                                                                       108812, 66252, 109487,
                                                                       39102, 39432, 71002,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 116762, 0, 3,
                                                                       110162, 67152, 110987,
                                                                       40092, 40488, 71552,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 117752, 0, 3,
                                                                       110987, 67702, 111812,
                                                                       40488, 40884, 72212,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 118742, 0, 3,
                                                                       111812, 68252, 112637,
                                                                       40884, 41280, 72872,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 119732, 0, 3,
                                                                       113462, 69352, 114287,
                                                                       42072, 42468, 73532,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 120722, 0, 3,
                                                                       114287, 69902, 115112,
                                                                       42468, 42864, 74192,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 121712, 0, 3,
                                                                       115112, 70452, 115937,
                                                                       42864, 43260, 74852,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122702, 3, 44052,
                                                                       44062, 75542, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122723, 3, 44062,
                                                                       44072, 75557, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122744, 3, 44072,
                                                                       44082, 75572, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122765, 3, 44082,
                                                                       44092, 75587, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122786, 3, 44092,
                                                                       44102, 75602, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122807, 3, 44102,
                                                                       44112, 75617, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122828, 3, 44112,
                                                                       44122, 75632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122849, 3, 44122,
                                                                       44132, 75647, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122870, 3, 44132,
                                                                       44142, 75662, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122891, 3, 44142,
                                                                       44152, 75677, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122912, 3, 44152,
                                                                       44162, 75692, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122933, 3, 44182,
                                                                       44192, 75737, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122954, 3, 44192,
                                                                       44202, 75752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122975, 3, 44202,
                                                                       44212, 75767, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 122996, 3, 44212,
                                                                       44222, 75782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 123017, 3, 44222,
                                                                       44232, 75797, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 123038, 3, 44232,
                                                                       44242, 75812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 123059, 3, 44242,
                                                                       44252, 75827, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 123080, 3, 44252,
                                                                       44262, 75842, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 123101, 3, 44262,
                                                                       44272, 75857, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 123122, 3, 44272,
                                                                       44282, 75872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 123143, 3, 44282,
                                                                       44292, 75887, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123164, 0, 3,
                                                                       122702, 75542, 122723,
                                                                       44312, 44342, 75992,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123227, 0, 3,
                                                                       122723, 75557, 122744,
                                                                       44342, 44372, 76037,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123290, 0, 3,
                                                                       122744, 75572, 122765,
                                                                       44372, 44402, 76082,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123353, 0, 3,
                                                                       122765, 75587, 122786,
                                                                       44402, 44432, 76127,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123416, 0, 3,
                                                                       122786, 75602, 122807,
                                                                       44432, 44462, 76172,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123479, 0, 3,
                                                                       122807, 75617, 122828,
                                                                       44462, 44492, 76217,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123542, 0, 3,
                                                                       122828, 75632, 122849,
                                                                       44492, 44522, 76262,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123605, 0, 3,
                                                                       122849, 75647, 122870,
                                                                       44522, 44552, 76307,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123668, 0, 3,
                                                                       122870, 75662, 122891,
                                                                       44552, 44582, 76352,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123731, 0, 3,
                                                                       122891, 75677, 122912,
                                                                       44582, 44612, 76397,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123794, 0, 3,
                                                                       122933, 75737, 122954,
                                                                       44672, 44702, 76532,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123857, 0, 3,
                                                                       122954, 75752, 122975,
                                                                       44702, 44732, 76577,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123920, 0, 3,
                                                                       122975, 75767, 122996,
                                                                       44732, 44762, 76622,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123983, 0, 3,
                                                                       122996, 75782, 123017,
                                                                       44762, 44792, 76667,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 124046, 0, 3,
                                                                       123017, 75797, 123038,
                                                                       44792, 44822, 76712,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 124109, 0, 3,
                                                                       123038, 75812, 123059,
                                                                       44822, 44852, 76757,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 124172, 0, 3,
                                                                       123059, 75827, 123080,
                                                                       44852, 44882, 76802,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 124235, 0, 3,
                                                                       123080, 75842, 123101,
                                                                       44882, 44912, 76847,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 124298, 0, 3,
                                                                       123101, 75857, 123122,
                                                                       44912, 44942, 76892,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 124361, 0, 3,
                                                                       123122, 75872, 123143,
                                                                       44942, 44972, 76937,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124424, 0, 3,
                                                                       123164, 75992, 123227,
                                                                       45032, 45092, 77162,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124550, 0, 3,
                                                                       123227, 76037, 123290,
                                                                       45092, 45152, 77252,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124676, 0, 3,
                                                                       123290, 76082, 123353,
                                                                       45152, 45212, 77342,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124802, 0, 3,
                                                                       123353, 76127, 123416,
                                                                       45212, 45272, 77432,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124928, 0, 3,
                                                                       123416, 76172, 123479,
                                                                       45272, 45332, 77522,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125054, 0, 3,
                                                                       123479, 76217, 123542,
                                                                       45332, 45392, 77612,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125180, 0, 3,
                                                                       123542, 76262, 123605,
                                                                       45392, 45452, 77702,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125306, 0, 3,
                                                                       123605, 76307, 123668,
                                                                       45452, 45512, 77792,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125432, 0, 3,
                                                                       123668, 76352, 123731,
                                                                       45512, 45572, 77882,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125558, 0, 3,
                                                                       123794, 76532, 123857,
                                                                       45692, 45752, 78152,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125684, 0, 3,
                                                                       123857, 76577, 123920,
                                                                       45752, 45812, 78242,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125810, 0, 3,
                                                                       123920, 76622, 123983,
                                                                       45812, 45872, 78332,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125936, 0, 3,
                                                                       123983, 76667, 124046,
                                                                       45872, 45932, 78422,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 126062, 0, 3,
                                                                       124046, 76712, 124109,
                                                                       45932, 45992, 78512,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 126188, 0, 3,
                                                                       124109, 76757, 124172,
                                                                       45992, 46052, 78602,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 126314, 0, 3,
                                                                       124172, 76802, 124235,
                                                                       46052, 46112, 78692,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 126440, 0, 3,
                                                                       124235, 76847, 124298,
                                                                       46112, 46172, 78782,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 126566, 0, 3,
                                                                       124298, 76892, 124361,
                                                                       46172, 46232, 78872,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 126692, 0, 3,
                                                                       124424, 77162, 124550,
                                                                       46352, 46452, 79262,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 126902, 0, 3,
                                                                       124550, 77252, 124676,
                                                                       46452, 46552, 79412,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 127112, 0, 3,
                                                                       124676, 77342, 124802,
                                                                       46552, 46652, 79562,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 127322, 0, 3,
                                                                       124802, 77432, 124928,
                                                                       46652, 46752, 79712,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 127532, 0, 3,
                                                                       124928, 77522, 125054,
                                                                       46752, 46852, 79862,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 127742, 0, 3,
                                                                       125054, 77612, 125180,
                                                                       46852, 46952, 80012,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 127952, 0, 3,
                                                                       125180, 77702, 125306,
                                                                       46952, 47052, 80162,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 128162, 0, 3,
                                                                       125306, 77792, 125432,
                                                                       47052, 47152, 80312,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 128372, 0, 3,
                                                                       125558, 78152, 125684,
                                                                       47352, 47452, 80762,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 128582, 0, 3,
                                                                       125684, 78242, 125810,
                                                                       47452, 47552, 80912,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 128792, 0, 3,
                                                                       125810, 78332, 125936,
                                                                       47552, 47652, 81062,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 129002, 0, 3,
                                                                       125936, 78422, 126062,
                                                                       47652, 47752, 81212,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 129212, 0, 3,
                                                                       126062, 78512, 126188,
                                                                       47752, 47852, 81362,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 129422, 0, 3,
                                                                       126188, 78602, 126314,
                                                                       47852, 47952, 81512,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 129632, 0, 3,
                                                                       126314, 78692, 126440,
                                                                       47952, 48052, 81662,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 129842, 0, 3,
                                                                       126440, 78782, 126566,
                                                                       48052, 48152, 81812,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 130052, 0, 3,
                                                                       126692, 79262, 126902,
                                                                       48352, 48502, 82412,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 130367, 0, 3,
                                                                       126902, 79412, 127112,
                                                                       48502, 48652, 82637,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 130682, 0, 3,
                                                                       127112, 79562, 127322,
                                                                       48652, 48802, 82862,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 130997, 0, 3,
                                                                       127322, 79712, 127532,
                                                                       48802, 48952, 83087,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 131312, 0, 3,
                                                                       127532, 79862, 127742,
                                                                       48952, 49102, 83312,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 131627, 0, 3,
                                                                       127742, 80012, 127952,
                                                                       49102, 49252, 83537,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 131942, 0, 3,
                                                                       127952, 80162, 128162,
                                                                       49252, 49402, 83762,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 132257, 0, 3,
                                                                       128372, 80762, 128582,
                                                                       49702, 49852, 84437,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 132572, 0, 3,
                                                                       128582, 80912, 128792,
                                                                       49852, 50002, 84662,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 132887, 0, 3,
                                                                       128792, 81062, 129002,
                                                                       50002, 50152, 84887,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 133202, 0, 3,
                                                                       129002, 81212, 129212,
                                                                       50152, 50302, 85112,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 133517, 0, 3,
                                                                       129212, 81362, 129422,
                                                                       50302, 50452, 85337,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 133832, 0, 3,
                                                                       129422, 81512, 129632,
                                                                       50452, 50602, 85562,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 134147, 0, 3,
                                                                       129632, 81662, 129842,
                                                                       50602, 50752, 85787,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 134462, 0, 3,
                                                                       130052, 82412, 130367,
                                                                       51052, 51262, 86642,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 134903, 0, 3,
                                                                       130367, 82637, 130682,
                                                                       51262, 51472, 86957,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 135344, 0, 3,
                                                                       130682, 82862, 130997,
                                                                       51472, 51682, 87272,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 135785, 0, 3,
                                                                       130997, 83087, 131312,
                                                                       51682, 51892, 87587,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 136226, 0, 3,
                                                                       131312, 83312, 131627,
                                                                       51892, 52102, 87902,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 136667, 0, 3,
                                                                       131627, 83537, 131942,
                                                                       52102, 52312, 88217,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 137108, 0, 3,
                                                                       132257, 84437, 132572,
                                                                       52732, 52942, 89162,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 137549, 0, 3,
                                                                       132572, 84662, 132887,
                                                                       52942, 53152, 89477,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 137990, 0, 3,
                                                                       132887, 84887, 133202,
                                                                       53152, 53362, 89792,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 138431, 0, 3,
                                                                       133202, 85112, 133517,
                                                                       53362, 53572, 90107,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 138872, 0, 3,
                                                                       133517, 85337, 133832,
                                                                       53572, 53782, 90422,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 139313, 0, 3,
                                                                       133832, 85562, 134147,
                                                                       53782, 53992, 90737,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 139754, 0, 3,
                                                                       134462, 86642, 134903,
                                                                       54412, 54692, 91892,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 140342, 0, 3,
                                                                       134903, 86957, 135344,
                                                                       54692, 54972, 92312,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 140930, 0, 3,
                                                                       135344, 87272, 135785,
                                                                       54972, 55252, 92732,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 141518, 0, 3,
                                                                       135785, 87587, 136226,
                                                                       55252, 55532, 93152,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 142106, 0, 3,
                                                                       136226, 87902, 136667,
                                                                       55532, 55812, 93572,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 142694, 0, 3,
                                                                       137108, 89162, 137549,
                                                                       56372, 56652, 94832,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 143282, 0, 3,
                                                                       137549, 89477, 137990,
                                                                       56652, 56932, 95252,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 143870, 0, 3,
                                                                       137990, 89792, 138431,
                                                                       56932, 57212, 95672,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 144458, 0, 3,
                                                                       138431, 90107, 138872,
                                                                       57212, 57492, 96092,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 145046, 0, 3,
                                                                       138872, 90422, 139313,
                                                                       57492, 57772, 96512,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 145634, 0, 3,
                                                                       139754, 91892, 140342,
                                                                       58332, 58692, 98012,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 146390, 0, 3,
                                                                       140342, 92312, 140930,
                                                                       58692, 59052, 98552,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 147146, 0, 3,
                                                                       140930, 92732, 141518,
                                                                       59052, 59412, 99092,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 147902, 0, 3,
                                                                       141518, 93152, 142106,
                                                                       59412, 59772, 99632,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 148658, 0, 3,
                                                                       142694, 94832, 143282,
                                                                       60492, 60852, 101252,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 149414, 0, 3,
                                                                       143282, 95252, 143870,
                                                                       60852, 61212, 101792,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 150170, 0, 3,
                                                                       143870, 95672, 144458,
                                                                       61212, 61572, 102332,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 150926, 0, 3,
                                                                       144458, 96092, 145046,
                                                                       61572, 61932, 102872,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 151682, 0, 3,
                                                                       145634, 98012, 146390,
                                                                       62652, 63102, 104762,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 152627, 0, 3,
                                                                       146390, 98552, 147146,
                                                                       63102, 63552, 105437,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 153572, 0, 3,
                                                                       147146, 99092, 147902,
                                                                       63552, 64002, 106112,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 154517, 0, 3,
                                                                       148658, 101252, 149414,
                                                                       64902, 65352, 108137,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 155462, 0, 3,
                                                                       149414, 101792, 150170,
                                                                       65352, 65802, 108812,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 156407, 0, 3,
                                                                       150170, 102332, 150926,
                                                                       65802, 66252, 109487,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 157352, 0, 3,
                                                                       151682, 104762, 152627,
                                                                       67152, 67702, 111812,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 158507, 0, 3,
                                                                       152627, 105437, 153572,
                                                                       67702, 68252, 112637,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 159662, 0, 3,
                                                                       154517, 108137, 155462,
                                                                       69352, 69902, 115112,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 160817, 0, 3,
                                                                       155462, 108812, 156407,
                                                                       69902, 70452, 115937,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 161972, 0, 3,
                                                                       157352, 111812, 158507,
                                                                       71552, 72212, 118742,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 163358, 0, 3,
                                                                       159662, 115112, 160817,
                                                                       73532, 74192, 121712,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 164744, 3, 75512,
                                                                       75527, 122702, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 164772, 3, 75527,
                                                                       75542, 122723, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 164800, 3, 75542,
                                                                       75557, 122744, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 164828, 3, 75557,
                                                                       75572, 122765, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 164856, 3, 75572,
                                                                       75587, 122786, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 164884, 3, 75587,
                                                                       75602, 122807, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 164912, 3, 75602,
                                                                       75617, 122828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 164940, 3, 75617,
                                                                       75632, 122849, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 164968, 3, 75632,
                                                                       75647, 122870, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 164996, 3, 75647,
                                                                       75662, 122891, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165024, 3, 75662,
                                                                       75677, 122912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165052, 3, 75707,
                                                                       75722, 122933, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165080, 3, 75722,
                                                                       75737, 122954, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165108, 3, 75737,
                                                                       75752, 122975, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165136, 3, 75752,
                                                                       75767, 122996, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165164, 3, 75767,
                                                                       75782, 123017, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165192, 3, 75782,
                                                                       75797, 123038, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165220, 3, 75797,
                                                                       75812, 123059, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165248, 3, 75812,
                                                                       75827, 123080, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165276, 3, 75827,
                                                                       75842, 123101, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165304, 3, 75842,
                                                                       75857, 123122, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 165332, 3, 75857,
                                                                       75872, 123143, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 165360, 0, 3,
                                                                       164744, 122702, 164772,
                                                                       75902, 75947, 123164,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 165444, 0, 3,
                                                                       164772, 122723, 164800,
                                                                       75947, 75992, 123227,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 165528, 0, 3,
                                                                       164800, 122744, 164828,
                                                                       75992, 76037, 123290,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 165612, 0, 3,
                                                                       164828, 122765, 164856,
                                                                       76037, 76082, 123353,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 165696, 0, 3,
                                                                       164856, 122786, 164884,
                                                                       76082, 76127, 123416,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 165780, 0, 3,
                                                                       164884, 122807, 164912,
                                                                       76127, 76172, 123479,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 165864, 0, 3,
                                                                       164912, 122828, 164940,
                                                                       76172, 76217, 123542,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 165948, 0, 3,
                                                                       164940, 122849, 164968,
                                                                       76217, 76262, 123605,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166032, 0, 3,
                                                                       164968, 122870, 164996,
                                                                       76262, 76307, 123668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166116, 0, 3,
                                                                       164996, 122891, 165024,
                                                                       76307, 76352, 123731,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166200, 0, 3,
                                                                       165052, 122933, 165080,
                                                                       76442, 76487, 123794,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166284, 0, 3,
                                                                       165080, 122954, 165108,
                                                                       76487, 76532, 123857,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166368, 0, 3,
                                                                       165108, 122975, 165136,
                                                                       76532, 76577, 123920,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166452, 0, 3,
                                                                       165136, 122996, 165164,
                                                                       76577, 76622, 123983,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166536, 0, 3,
                                                                       165164, 123017, 165192,
                                                                       76622, 76667, 124046,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166620, 0, 3,
                                                                       165192, 123038, 165220,
                                                                       76667, 76712, 124109,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166704, 0, 3,
                                                                       165220, 123059, 165248,
                                                                       76712, 76757, 124172,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166788, 0, 3,
                                                                       165248, 123080, 165276,
                                                                       76757, 76802, 124235,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166872, 0, 3,
                                                                       165276, 123101, 165304,
                                                                       76802, 76847, 124298,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 166956, 0, 3,
                                                                       165304, 123122, 165332,
                                                                       76847, 76892, 124361,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 167040, 0, 3,
                                                                       165360, 123164, 165444,
                                                                       76982, 77072, 124424,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 167208, 0, 3,
                                                                       165444, 123227, 165528,
                                                                       77072, 77162, 124550,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 167376, 0, 3,
                                                                       165528, 123290, 165612,
                                                                       77162, 77252, 124676,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 167544, 0, 3,
                                                                       165612, 123353, 165696,
                                                                       77252, 77342, 124802,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 167712, 0, 3,
                                                                       165696, 123416, 165780,
                                                                       77342, 77432, 124928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 167880, 0, 3,
                                                                       165780, 123479, 165864,
                                                                       77432, 77522, 125054,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 168048, 0, 3,
                                                                       165864, 123542, 165948,
                                                                       77522, 77612, 125180,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 168216, 0, 3,
                                                                       165948, 123605, 166032,
                                                                       77612, 77702, 125306,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 168384, 0, 3,
                                                                       166032, 123668, 166116,
                                                                       77702, 77792, 125432,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 168552, 0, 3,
                                                                       166200, 123794, 166284,
                                                                       77972, 78062, 125558,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 168720, 0, 3,
                                                                       166284, 123857, 166368,
                                                                       78062, 78152, 125684,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 168888, 0, 3,
                                                                       166368, 123920, 166452,
                                                                       78152, 78242, 125810,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 169056, 0, 3,
                                                                       166452, 123983, 166536,
                                                                       78242, 78332, 125936,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 169224, 0, 3,
                                                                       166536, 124046, 166620,
                                                                       78332, 78422, 126062,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 169392, 0, 3,
                                                                       166620, 124109, 166704,
                                                                       78422, 78512, 126188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 169560, 0, 3,
                                                                       166704, 124172, 166788,
                                                                       78512, 78602, 126314,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 169728, 0, 3,
                                                                       166788, 124235, 166872,
                                                                       78602, 78692, 126440,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 169896, 0, 3,
                                                                       166872, 124298, 166956,
                                                                       78692, 78782, 126566,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 170064, 0, 3,
                                                                       167040, 124424, 167208,
                                                                       78962, 79112, 126692,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 170344, 0, 3,
                                                                       167208, 124550, 167376,
                                                                       79112, 79262, 126902,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 170624, 0, 3,
                                                                       167376, 124676, 167544,
                                                                       79262, 79412, 127112,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 170904, 0, 3,
                                                                       167544, 124802, 167712,
                                                                       79412, 79562, 127322,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 171184, 0, 3,
                                                                       167712, 124928, 167880,
                                                                       79562, 79712, 127532,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 171464, 0, 3,
                                                                       167880, 125054, 168048,
                                                                       79712, 79862, 127742,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 171744, 0, 3,
                                                                       168048, 125180, 168216,
                                                                       79862, 80012, 127952,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 172024, 0, 3,
                                                                       168216, 125306, 168384,
                                                                       80012, 80162, 128162,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 172304, 0, 3,
                                                                       168552, 125558, 168720,
                                                                       80462, 80612, 128372,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 172584, 0, 3,
                                                                       168720, 125684, 168888,
                                                                       80612, 80762, 128582,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 172864, 0, 3,
                                                                       168888, 125810, 169056,
                                                                       80762, 80912, 128792,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 173144, 0, 3,
                                                                       169056, 125936, 169224,
                                                                       80912, 81062, 129002,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 173424, 0, 3,
                                                                       169224, 126062, 169392,
                                                                       81062, 81212, 129212,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 173704, 0, 3,
                                                                       169392, 126188, 169560,
                                                                       81212, 81362, 129422,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 173984, 0, 3,
                                                                       169560, 126314, 169728,
                                                                       81362, 81512, 129632,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 174264, 0, 3,
                                                                       169728, 126440, 169896,
                                                                       81512, 81662, 129842,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 174544, 0, 3,
                                                                       170064, 126692, 170344,
                                                                       81962, 82187, 130052,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 174964, 0, 3,
                                                                       170344, 126902, 170624,
                                                                       82187, 82412, 130367,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 175384, 0, 3,
                                                                       170624, 127112, 170904,
                                                                       82412, 82637, 130682,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 175804, 0, 3,
                                                                       170904, 127322, 171184,
                                                                       82637, 82862, 130997,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 176224, 0, 3,
                                                                       171184, 127532, 171464,
                                                                       82862, 83087, 131312,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 176644, 0, 3,
                                                                       171464, 127742, 171744,
                                                                       83087, 83312, 131627,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 177064, 0, 3,
                                                                       171744, 127952, 172024,
                                                                       83312, 83537, 131942,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 177484, 0, 3,
                                                                       172304, 128372, 172584,
                                                                       83987, 84212, 132257,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 177904, 0, 3,
                                                                       172584, 128582, 172864,
                                                                       84212, 84437, 132572,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 178324, 0, 3,
                                                                       172864, 128792, 173144,
                                                                       84437, 84662, 132887,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 178744, 0, 3,
                                                                       173144, 129002, 173424,
                                                                       84662, 84887, 133202,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 179164, 0, 3,
                                                                       173424, 129212, 173704,
                                                                       84887, 85112, 133517,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 179584, 0, 3,
                                                                       173704, 129422, 173984,
                                                                       85112, 85337, 133832,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 180004, 0, 3,
                                                                       173984, 129632, 174264,
                                                                       85337, 85562, 134147,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 180424, 0, 3,
                                                                       174544, 130052, 174964,
                                                                       86012, 86327, 134462,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 181012, 0, 3,
                                                                       174964, 130367, 175384,
                                                                       86327, 86642, 134903,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 181600, 0, 3,
                                                                       175384, 130682, 175804,
                                                                       86642, 86957, 135344,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 182188, 0, 3,
                                                                       175804, 130997, 176224,
                                                                       86957, 87272, 135785,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 182776, 0, 3,
                                                                       176224, 131312, 176644,
                                                                       87272, 87587, 136226,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 183364, 0, 3,
                                                                       176644, 131627, 177064,
                                                                       87587, 87902, 136667,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 183952, 0, 3,
                                                                       177484, 132257, 177904,
                                                                       88532, 88847, 137108,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 184540, 0, 3,
                                                                       177904, 132572, 178324,
                                                                       88847, 89162, 137549,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 185128, 0, 3,
                                                                       178324, 132887, 178744,
                                                                       89162, 89477, 137990,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 185716, 0, 3,
                                                                       178744, 133202, 179164,
                                                                       89477, 89792, 138431,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 186304, 0, 3,
                                                                       179164, 133517, 179584,
                                                                       89792, 90107, 138872,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 186892, 0, 3,
                                                                       179584, 133832, 180004,
                                                                       90107, 90422, 139313,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 187480, 0, 3,
                                                                       180424, 134462, 181012,
                                                                       91052, 91472, 139754,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 188264, 0, 3,
                                                                       181012, 134903, 181600,
                                                                       91472, 91892, 140342,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 189048, 0, 3,
                                                                       181600, 135344, 182188,
                                                                       91892, 92312, 140930,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 189832, 0, 3,
                                                                       182188, 135785, 182776,
                                                                       92312, 92732, 141518,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 190616, 0, 3,
                                                                       182776, 136226, 183364,
                                                                       92732, 93152, 142106,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 191400, 0, 3,
                                                                       183952, 137108, 184540,
                                                                       93992, 94412, 142694,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 192184, 0, 3,
                                                                       184540, 137549, 185128,
                                                                       94412, 94832, 143282,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 192968, 0, 3,
                                                                       185128, 137990, 185716,
                                                                       94832, 95252, 143870,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 193752, 0, 3,
                                                                       185716, 138431, 186304,
                                                                       95252, 95672, 144458,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 194536, 0, 3,
                                                                       186304, 138872, 186892,
                                                                       95672, 96092, 145046,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 195320, 0, 3,
                                                                       187480, 139754, 188264,
                                                                       96932, 97472, 145634,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 196328, 0, 3,
                                                                       188264, 140342, 189048,
                                                                       97472, 98012, 146390,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 197336, 0, 3,
                                                                       189048, 140930, 189832,
                                                                       98012, 98552, 147146,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 198344, 0, 3,
                                                                       189832, 141518, 190616,
                                                                       98552, 99092, 147902,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 199352, 0, 3,
                                                                       191400, 142694, 192184,
                                                                       100172, 100712, 148658,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 200360, 0, 3,
                                                                       192184, 143282, 192968,
                                                                       100712, 101252, 149414,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 201368, 0, 3,
                                                                       192968, 143870, 193752,
                                                                       101252, 101792, 150170,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 202376, 0, 3,
                                                                       193752, 144458, 194536,
                                                                       101792, 102332, 150926,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 203384, 0, 3,
                                                                       195320, 145634, 196328,
                                                                       103412, 104087, 151682,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 204644, 0, 3,
                                                                       196328, 146390, 197336,
                                                                       104087, 104762, 152627,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 205904, 0, 3,
                                                                       197336, 147146, 198344,
                                                                       104762, 105437, 153572,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 207164, 0, 3,
                                                                       199352, 148658, 200360,
                                                                       106787, 107462, 154517,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 208424, 0, 3,
                                                                       200360, 149414, 201368,
                                                                       107462, 108137, 155462,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 209684, 0, 3,
                                                                       201368, 150170, 202376,
                                                                       108137, 108812, 156407,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 210944, 0, 3,
                                                                       203384, 151682, 204644,
                                                                       110162, 110987, 157352,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 212484, 0, 3,
                                                                       204644, 152627, 205904,
                                                                       110987, 111812, 158507,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 214024, 0, 3,
                                                                       207164, 154517, 208424,
                                                                       113462, 114287, 159662,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 215564, 0, 3,
                                                                       208424, 155462, 209684,
                                                                       114287, 115112, 160817,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 217104, 0, 3,
                                                                       210944, 157352, 212484,
                                                                       116762, 117752, 161972,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 218952, 0, 3,
                                                                       214024, 159662, 215564,
                                                                       119732, 120722, 163358,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 220800, 180424, 588, ncols);

                    simdfunc::contract_primitives(buffer, 221661, 183952, 588, ncols);

                    simdfunc::contract_primitives(buffer, 222522, 187480, 784, ncols);

                    simdfunc::contract_primitives(buffer, 223670, 191400, 784, ncols);

                    simdfunc::contract_primitives(buffer, 224818, 195320, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 226294, 199352, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 227770, 203384, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 229615, 207164, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 231460, 210944, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 233715, 214024, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 235970, 217104, 1848, ncols);

                    simdfunc::contract_primitives(buffer, 238676, 218952, 1848, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 221388, 220800, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 222249, 221661, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 223306, 222522, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 224454, 223670, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 225826, 224818, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 227302, 226294, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 229030, 227770, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 230875, 229615, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 233000, 231460, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 235255, 233715, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 237818, 235970, 66, 1, nmax);

        simdtrf::transform_i_inner(buffer, 240524, 238676, 66, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 241382, 221388, 223306, 13,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 242201, 222249, 224454, 13,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 243020, 223306, 225826, 13,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 244112, 224454, 227302, 13,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 245204, 225826, 229030, 13,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 246608, 227302, 230875, 13,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 248012, 229030, 233000, 13,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 249767, 230875, 235255, 13,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 251522, 233000, 237818, 13,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 253667, 235255, 240524, 13,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 255812, 241382, 243020, 13,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 257450, 242201, 244112, 13,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 259088, 243020, 245204, 13,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 261272, 244112, 246608, 13,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 263456, 245204, 248012, 13,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 266264, 246608, 249767, 13,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 269072, 248012, 251522, 13,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 272582, 249767, 253667, 13,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 276092, 255812, 259088, 13,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 278822, 257450, 261272, 13,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 281552, 259088, 263456, 13,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 285192, 261272, 266264, 13,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 288832, 263456, 269072, 13,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 293512, 266264, 272582, 13,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 298192, 276092, 281552, 13,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 302287, 278822, 285192, 13,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 306382, 281552, 288832, 13,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 311842, 285192, 293512, 13,
                                             nmax);

        simdtrf::compute_hrr_hh_out_of_first(buffer, coordinates, 317302, 298192, 306382, 13,
                                             nmax);

        simdtrf::compute_hrr_hh_out_of_first(buffer, coordinates, 323035, 302287, 311842, 13,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 328768, 323035, 21, 13, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 328768, 143, nmax);

        simdtrf::transform_h_inner(buffer, 328768, 317302, 21, 13, nmax);

        simdtrf::transform_h_outer(values + 1573 * nvalues + n * npairs, nvalues, buffer, 328768,
                                   143, nmax);
    }

    for (size_t m = 0; m < 3146; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
