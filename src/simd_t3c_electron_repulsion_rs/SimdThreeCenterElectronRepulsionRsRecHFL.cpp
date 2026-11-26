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


#include "SimdThreeCenterElectronRepulsionRsRecHFL.hpp"

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
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_hfl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hfl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 293583, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2618 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 293583, 249158, 15355, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3638, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3641, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3644, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3647, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3650, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3653, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3656, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3659, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3662, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3665, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3668, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3671, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3674, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3677, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3680, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3683, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3686, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3689, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3692, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3695, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3698, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3701, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3704, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3707, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3710, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3713, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3716, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3719, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3722, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3725, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3728, 3, 9, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3737, 3, 10, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3746, 3, 11, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3755, 3, 12, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3764, 3, 13, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3773, 3, 14, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3782, 3, 15, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3791, 3, 16, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3800, 3, 17, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3809, 3, 18, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3818, 3, 19, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3827, 3, 20, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3836, 3, 21, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3845, 3, 22, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3854, 3, 27, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3863, 3, 28, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3872, 3, 29, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3881, 3, 30, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3890, 3, 31, 108,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3899, 3, 32, 111,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3908, 3, 33, 114,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3917, 3, 34, 117,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3926, 3, 35, 120,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3935, 3, 36, 123,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3944, 3, 37, 126,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3953, 3, 38, 129,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3962, 3, 39, 132,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3971, 3, 40, 135,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3980, 3, 48, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3998, 3, 51, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4016, 3, 54, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4034, 3, 57, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4052, 3, 60, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4070, 3, 63, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4088, 3, 66, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4106, 3, 69, 192,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4124, 3, 72, 198,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4142, 3, 75, 204,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4160, 3, 78, 210,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4178, 3, 81, 216,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4196, 3, 84, 222,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4214, 3, 96, 240,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4232, 3, 99, 246,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4250, 3, 102, 252,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4268, 3, 105, 258,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4286, 3, 108, 264,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4304, 3, 111, 270,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4322, 3, 114, 276,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4340, 3, 117, 282,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4358, 3, 120, 288,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4376, 3, 123, 294,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4394, 3, 126, 300,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4412, 3, 129, 306,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4430, 3, 132, 312,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4448, 3, 150, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4478, 3, 156, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4508, 3, 162, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4538, 3, 168, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4568, 3, 174, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4598, 3, 180, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4628, 3, 186, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4658, 3, 192, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4688, 3, 198, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4718, 3, 204, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4748, 3, 210, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4778, 3, 216, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4808, 3, 240, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4838, 3, 246, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4868, 3, 252, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4898, 3, 258, 508,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4928, 3, 264, 518,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4958, 3, 270, 528,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4988, 3, 276, 538,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5018, 3, 282, 548,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5048, 3, 288, 558,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5078, 3, 294, 568,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5108, 3, 300, 578,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5138, 3, 306, 588,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5168, 3, 338, 628,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5213, 3, 348, 643,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5258, 3, 358, 658,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5303, 3, 368, 673,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5348, 3, 378, 688,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5393, 3, 388, 703,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5438, 3, 398, 718,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5483, 3, 408, 733,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5528, 3, 418, 748,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5573, 3, 428, 763,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5618, 3, 438, 778,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5663, 3, 478, 823,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5708, 3, 488, 838,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5753, 3, 498, 853,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5798, 3, 508, 868,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5843, 3, 518, 883,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5888, 3, 528, 898,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5933, 3, 538, 913,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5978, 3, 548, 928,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6023, 3, 558, 943,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6068, 3, 568, 958,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6113, 3, 578, 973,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6158, 3, 628,
                                                                       1030, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6221, 3, 643,
                                                                       1051, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6284, 3, 658,
                                                                       1072, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6347, 3, 673,
                                                                       1093, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6410, 3, 688,
                                                                       1114, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6473, 3, 703,
                                                                       1135, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6536, 3, 718,
                                                                       1156, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6599, 3, 733,
                                                                       1177, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6662, 3, 748,
                                                                       1198, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6725, 3, 763,
                                                                       1219, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6788, 3, 823,
                                                                       1282, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6851, 3, 838,
                                                                       1303, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6914, 3, 853,
                                                                       1324, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6977, 3, 868,
                                                                       1345, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7040, 3, 883,
                                                                       1366, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7103, 3, 898,
                                                                       1387, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7166, 3, 913,
                                                                       1408, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7229, 3, 928,
                                                                       1429, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7292, 3, 943,
                                                                       1450, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7355, 3, 958,
                                                                       1471, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7418, 3, 1030,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7502, 3, 1051,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7586, 3, 1072,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7670, 3, 1093,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7754, 3, 1114,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7838, 3, 1135,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7922, 3, 1156,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8006, 3, 1177,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8090, 3, 1198,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8174, 3, 1282,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8258, 3, 1303,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8342, 3, 1324,
                                                                       1912, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8426, 3, 1345,
                                                                       1940, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8510, 3, 1366,
                                                                       1968, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8594, 3, 1387,
                                                                       1996, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8678, 3, 1408,
                                                                       2024, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8762, 3, 1429,
                                                                       2052, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8846, 3, 1450,
                                                                       2080, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8930, 3, 1548,
                                                                       2180, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9038, 3, 1576,
                                                                       2216, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9146, 3, 1604,
                                                                       2252, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9254, 3, 1632,
                                                                       2288, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9362, 3, 1660,
                                                                       2324, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9470, 3, 1688,
                                                                       2360, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9578, 3, 1716,
                                                                       2396, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9686, 3, 1744,
                                                                       2432, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9794, 3, 1856,
                                                                       2540, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9902, 3, 1884,
                                                                       2576, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10010, 3, 1912,
                                                                       2612, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10118, 3, 1940,
                                                                       2648, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10226, 3, 1968,
                                                                       2684, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10334, 3, 1996,
                                                                       2720, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10442, 3, 2024,
                                                                       2756, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10550, 3, 2052,
                                                                       2792, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10658, 3, 2180,
                                                                       2918, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10793, 3, 2216,
                                                                       2963, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10928, 3, 2252,
                                                                       3008, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11063, 3, 2288,
                                                                       3053, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11198, 3, 2324,
                                                                       3098, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11333, 3, 2360,
                                                                       3143, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11468, 3, 2396,
                                                                       3188, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11603, 3, 2540,
                                                                       3323, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11738, 3, 2576,
                                                                       3368, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11873, 3, 2612,
                                                                       3413, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12008, 3, 2648,
                                                                       3458, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12143, 3, 2684,
                                                                       3503, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12278, 3, 2720,
                                                                       3548, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12413, 3, 2756,
                                                                       3593, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12548, 3, 7, 8,
                                                                       3638, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12554, 3, 8, 9,
                                                                       3641, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12560, 3, 9, 10,
                                                                       3644, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12566, 3, 10, 11,
                                                                       3647, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12572, 3, 11, 12,
                                                                       3650, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12578, 3, 12, 13,
                                                                       3653, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12584, 3, 13, 14,
                                                                       3656, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12590, 3, 14, 15,
                                                                       3659, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12596, 3, 15, 16,
                                                                       3662, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12602, 3, 16, 17,
                                                                       3665, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12608, 3, 17, 18,
                                                                       3668, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12614, 3, 18, 19,
                                                                       3671, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12620, 3, 19, 20,
                                                                       3674, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12626, 3, 20, 21,
                                                                       3677, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12632, 3, 21, 22,
                                                                       3680, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12638, 3, 25, 26,
                                                                       3683, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12644, 3, 26, 27,
                                                                       3686, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12650, 3, 27, 28,
                                                                       3689, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12656, 3, 28, 29,
                                                                       3692, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12662, 3, 29, 30,
                                                                       3695, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12668, 3, 30, 31,
                                                                       3698, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12674, 3, 31, 32,
                                                                       3701, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12680, 3, 32, 33,
                                                                       3704, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12686, 3, 33, 34,
                                                                       3707, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12692, 3, 34, 35,
                                                                       3710, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12698, 3, 35, 36,
                                                                       3713, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12704, 3, 36, 37,
                                                                       3716, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12710, 3, 37, 38,
                                                                       3719, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12716, 3, 38, 39,
                                                                       3722, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12722, 3, 39, 40,
                                                                       3725, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12728, 0, 3,
                                                                       12548, 3638, 12554, 3728,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12746, 0, 3,
                                                                       12554, 3641, 12560, 3737,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12764, 0, 3,
                                                                       12560, 3644, 12566, 3746,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12782, 0, 3,
                                                                       12566, 3647, 12572, 3755,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12800, 0, 3,
                                                                       12572, 3650, 12578, 3764,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12818, 0, 3,
                                                                       12578, 3653, 12584, 3773,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12836, 0, 3,
                                                                       12584, 3656, 12590, 3782,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12854, 0, 3,
                                                                       12590, 3659, 12596, 3791,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12872, 0, 3,
                                                                       12596, 3662, 12602, 3800,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12890, 0, 3,
                                                                       12602, 3665, 12608, 3809,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12908, 0, 3,
                                                                       12608, 3668, 12614, 3818,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12926, 0, 3,
                                                                       12614, 3671, 12620, 3827,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12944, 0, 3,
                                                                       12620, 3674, 12626, 3836,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12962, 0, 3,
                                                                       12626, 3677, 12632, 3845,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12980, 0, 3,
                                                                       12638, 3683, 12644, 3854,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12998, 0, 3,
                                                                       12644, 3686, 12650, 3863,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13016, 0, 3,
                                                                       12650, 3689, 12656, 3872,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13034, 0, 3,
                                                                       12656, 3692, 12662, 3881,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13052, 0, 3,
                                                                       12662, 3695, 12668, 3890,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13070, 0, 3,
                                                                       12668, 3698, 12674, 3899,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13088, 0, 3,
                                                                       12674, 3701, 12680, 3908,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13106, 0, 3,
                                                                       12680, 3704, 12686, 3917,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13124, 0, 3,
                                                                       12686, 3707, 12692, 3926,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13142, 0, 3,
                                                                       12692, 3710, 12698, 3935,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13160, 0, 3,
                                                                       12698, 3713, 12704, 3944,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13178, 0, 3,
                                                                       12704, 3716, 12710, 3953,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13196, 0, 3,
                                                                       12710, 3719, 12716, 3962,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13214, 0, 3,
                                                                       12716, 3722, 12722, 3971,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13232, 0, 3,
                                                                       12728, 3728, 12746, 138,
                                                                       144, 3980, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13268, 0, 3,
                                                                       12746, 3737, 12764, 144,
                                                                       150, 3998, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13304, 0, 3,
                                                                       12764, 3746, 12782, 150,
                                                                       156, 4016, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13340, 0, 3,
                                                                       12782, 3755, 12800, 156,
                                                                       162, 4034, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13376, 0, 3,
                                                                       12800, 3764, 12818, 162,
                                                                       168, 4052, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13412, 0, 3,
                                                                       12818, 3773, 12836, 168,
                                                                       174, 4070, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13448, 0, 3,
                                                                       12836, 3782, 12854, 174,
                                                                       180, 4088, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13484, 0, 3,
                                                                       12854, 3791, 12872, 180,
                                                                       186, 4106, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13520, 0, 3,
                                                                       12872, 3800, 12890, 186,
                                                                       192, 4124, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13556, 0, 3,
                                                                       12890, 3809, 12908, 192,
                                                                       198, 4142, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13592, 0, 3,
                                                                       12908, 3818, 12926, 198,
                                                                       204, 4160, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13628, 0, 3,
                                                                       12926, 3827, 12944, 204,
                                                                       210, 4178, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13664, 0, 3,
                                                                       12944, 3836, 12962, 210,
                                                                       216, 4196, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13700, 0, 3,
                                                                       12980, 3854, 12998, 228,
                                                                       234, 4214, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13736, 0, 3,
                                                                       12998, 3863, 13016, 234,
                                                                       240, 4232, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13772, 0, 3,
                                                                       13016, 3872, 13034, 240,
                                                                       246, 4250, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13808, 0, 3,
                                                                       13034, 3881, 13052, 246,
                                                                       252, 4268, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13844, 0, 3,
                                                                       13052, 3890, 13070, 252,
                                                                       258, 4286, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13880, 0, 3,
                                                                       13070, 3899, 13088, 258,
                                                                       264, 4304, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13916, 0, 3,
                                                                       13088, 3908, 13106, 264,
                                                                       270, 4322, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13952, 0, 3,
                                                                       13106, 3917, 13124, 270,
                                                                       276, 4340, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13988, 0, 3,
                                                                       13124, 3926, 13142, 276,
                                                                       282, 4358, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14024, 0, 3,
                                                                       13142, 3935, 13160, 282,
                                                                       288, 4376, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14060, 0, 3,
                                                                       13160, 3944, 13178, 288,
                                                                       294, 4394, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14096, 0, 3,
                                                                       13178, 3953, 13196, 294,
                                                                       300, 4412, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14132, 0, 3,
                                                                       13196, 3962, 13214, 300,
                                                                       306, 4430, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14168, 0, 3,
                                                                       13232, 3980, 13268, 318,
                                                                       328, 4448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14228, 0, 3,
                                                                       13268, 3998, 13304, 328,
                                                                       338, 4478, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14288, 0, 3,
                                                                       13304, 4016, 13340, 338,
                                                                       348, 4508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14348, 0, 3,
                                                                       13340, 4034, 13376, 348,
                                                                       358, 4538, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14408, 0, 3,
                                                                       13376, 4052, 13412, 358,
                                                                       368, 4568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14468, 0, 3,
                                                                       13412, 4070, 13448, 368,
                                                                       378, 4598, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14528, 0, 3,
                                                                       13448, 4088, 13484, 378,
                                                                       388, 4628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14588, 0, 3,
                                                                       13484, 4106, 13520, 388,
                                                                       398, 4658, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14648, 0, 3,
                                                                       13520, 4124, 13556, 398,
                                                                       408, 4688, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14708, 0, 3,
                                                                       13556, 4142, 13592, 408,
                                                                       418, 4718, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14768, 0, 3,
                                                                       13592, 4160, 13628, 418,
                                                                       428, 4748, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14828, 0, 3,
                                                                       13628, 4178, 13664, 428,
                                                                       438, 4778, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14888, 0, 3,
                                                                       13700, 4214, 13736, 458,
                                                                       468, 4808, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14948, 0, 3,
                                                                       13736, 4232, 13772, 468,
                                                                       478, 4838, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15008, 0, 3,
                                                                       13772, 4250, 13808, 478,
                                                                       488, 4868, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15068, 0, 3,
                                                                       13808, 4268, 13844, 488,
                                                                       498, 4898, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15128, 0, 3,
                                                                       13844, 4286, 13880, 498,
                                                                       508, 4928, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15188, 0, 3,
                                                                       13880, 4304, 13916, 508,
                                                                       518, 4958, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15248, 0, 3,
                                                                       13916, 4322, 13952, 518,
                                                                       528, 4988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15308, 0, 3,
                                                                       13952, 4340, 13988, 528,
                                                                       538, 5018, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15368, 0, 3,
                                                                       13988, 4358, 14024, 538,
                                                                       548, 5048, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15428, 0, 3,
                                                                       14024, 4376, 14060, 548,
                                                                       558, 5078, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15488, 0, 3,
                                                                       14060, 4394, 14096, 558,
                                                                       568, 5108, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15548, 0, 3,
                                                                       14096, 4412, 14132, 568,
                                                                       578, 5138, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15608, 0, 3,
                                                                       14168, 4448, 14228, 598,
                                                                       613, 5168, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15698, 0, 3,
                                                                       14228, 4478, 14288, 613,
                                                                       628, 5213, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15788, 0, 3,
                                                                       14288, 4508, 14348, 628,
                                                                       643, 5258, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15878, 0, 3,
                                                                       14348, 4538, 14408, 643,
                                                                       658, 5303, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15968, 0, 3,
                                                                       14408, 4568, 14468, 658,
                                                                       673, 5348, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16058, 0, 3,
                                                                       14468, 4598, 14528, 673,
                                                                       688, 5393, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16148, 0, 3,
                                                                       14528, 4628, 14588, 688,
                                                                       703, 5438, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16238, 0, 3,
                                                                       14588, 4658, 14648, 703,
                                                                       718, 5483, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16328, 0, 3,
                                                                       14648, 4688, 14708, 718,
                                                                       733, 5528, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16418, 0, 3,
                                                                       14708, 4718, 14768, 733,
                                                                       748, 5573, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16508, 0, 3,
                                                                       14768, 4748, 14828, 748,
                                                                       763, 5618, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16598, 0, 3,
                                                                       14888, 4808, 14948, 793,
                                                                       808, 5663, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16688, 0, 3,
                                                                       14948, 4838, 15008, 808,
                                                                       823, 5708, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16778, 0, 3,
                                                                       15008, 4868, 15068, 823,
                                                                       838, 5753, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16868, 0, 3,
                                                                       15068, 4898, 15128, 838,
                                                                       853, 5798, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16958, 0, 3,
                                                                       15128, 4928, 15188, 853,
                                                                       868, 5843, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17048, 0, 3,
                                                                       15188, 4958, 15248, 868,
                                                                       883, 5888, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17138, 0, 3,
                                                                       15248, 4988, 15308, 883,
                                                                       898, 5933, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17228, 0, 3,
                                                                       15308, 5018, 15368, 898,
                                                                       913, 5978, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17318, 0, 3,
                                                                       15368, 5048, 15428, 913,
                                                                       928, 6023, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17408, 0, 3,
                                                                       15428, 5078, 15488, 928,
                                                                       943, 6068, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17498, 0, 3,
                                                                       15488, 5108, 15548, 943,
                                                                       958, 6113, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17588, 0, 3,
                                                                       15608, 5168, 15698, 988,
                                                                       1009, 6158, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17714, 0, 3,
                                                                       15698, 5213, 15788, 1009,
                                                                       1030, 6221, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17840, 0, 3,
                                                                       15788, 5258, 15878, 1030,
                                                                       1051, 6284, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17966, 0, 3,
                                                                       15878, 5303, 15968, 1051,
                                                                       1072, 6347, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18092, 0, 3,
                                                                       15968, 5348, 16058, 1072,
                                                                       1093, 6410, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18218, 0, 3,
                                                                       16058, 5393, 16148, 1093,
                                                                       1114, 6473, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18344, 0, 3,
                                                                       16148, 5438, 16238, 1114,
                                                                       1135, 6536, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18470, 0, 3,
                                                                       16238, 5483, 16328, 1135,
                                                                       1156, 6599, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18596, 0, 3,
                                                                       16328, 5528, 16418, 1156,
                                                                       1177, 6662, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18722, 0, 3,
                                                                       16418, 5573, 16508, 1177,
                                                                       1198, 6725, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18848, 0, 3,
                                                                       16598, 5663, 16688, 1240,
                                                                       1261, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18974, 0, 3,
                                                                       16688, 5708, 16778, 1261,
                                                                       1282, 6851, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19100, 0, 3,
                                                                       16778, 5753, 16868, 1282,
                                                                       1303, 6914, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19226, 0, 3,
                                                                       16868, 5798, 16958, 1303,
                                                                       1324, 6977, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19352, 0, 3,
                                                                       16958, 5843, 17048, 1324,
                                                                       1345, 7040, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19478, 0, 3,
                                                                       17048, 5888, 17138, 1345,
                                                                       1366, 7103, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19604, 0, 3,
                                                                       17138, 5933, 17228, 1366,
                                                                       1387, 7166, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19730, 0, 3,
                                                                       17228, 5978, 17318, 1387,
                                                                       1408, 7229, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19856, 0, 3,
                                                                       17318, 6023, 17408, 1408,
                                                                       1429, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19982, 0, 3,
                                                                       17408, 6068, 17498, 1429,
                                                                       1450, 7355, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20108, 0, 3,
                                                                       17588, 6158, 17714, 1492,
                                                                       1520, 7418, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20276, 0, 3,
                                                                       17714, 6221, 17840, 1520,
                                                                       1548, 7502, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20444, 0, 3,
                                                                       17840, 6284, 17966, 1548,
                                                                       1576, 7586, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20612, 0, 3,
                                                                       17966, 6347, 18092, 1576,
                                                                       1604, 7670, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20780, 0, 3,
                                                                       18092, 6410, 18218, 1604,
                                                                       1632, 7754, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20948, 0, 3,
                                                                       18218, 6473, 18344, 1632,
                                                                       1660, 7838, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21116, 0, 3,
                                                                       18344, 6536, 18470, 1660,
                                                                       1688, 7922, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21284, 0, 3,
                                                                       18470, 6599, 18596, 1688,
                                                                       1716, 8006, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21452, 0, 3,
                                                                       18596, 6662, 18722, 1716,
                                                                       1744, 8090, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21620, 0, 3,
                                                                       18848, 6788, 18974, 1800,
                                                                       1828, 8174, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21788, 0, 3,
                                                                       18974, 6851, 19100, 1828,
                                                                       1856, 8258, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21956, 0, 3,
                                                                       19100, 6914, 19226, 1856,
                                                                       1884, 8342, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22124, 0, 3,
                                                                       19226, 6977, 19352, 1884,
                                                                       1912, 8426, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22292, 0, 3,
                                                                       19352, 7040, 19478, 1912,
                                                                       1940, 8510, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22460, 0, 3,
                                                                       19478, 7103, 19604, 1940,
                                                                       1968, 8594, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22628, 0, 3,
                                                                       19604, 7166, 19730, 1968,
                                                                       1996, 8678, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22796, 0, 3,
                                                                       19730, 7229, 19856, 1996,
                                                                       2024, 8762, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22964, 0, 3,
                                                                       19856, 7292, 19982, 2024,
                                                                       2052, 8846, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23132, 0, 3,
                                                                       20108, 7418, 20276, 2108,
                                                                       2144, 8930, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23348, 0, 3,
                                                                       20276, 7502, 20444, 2144,
                                                                       2180, 9038, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23564, 0, 3,
                                                                       20444, 7586, 20612, 2180,
                                                                       2216, 9146, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23780, 0, 3,
                                                                       20612, 7670, 20780, 2216,
                                                                       2252, 9254, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23996, 0, 3,
                                                                       20780, 7754, 20948, 2252,
                                                                       2288, 9362, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24212, 0, 3,
                                                                       20948, 7838, 21116, 2288,
                                                                       2324, 9470, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24428, 0, 3,
                                                                       21116, 7922, 21284, 2324,
                                                                       2360, 9578, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24644, 0, 3,
                                                                       21284, 8006, 21452, 2360,
                                                                       2396, 9686, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24860, 0, 3,
                                                                       21620, 8174, 21788, 2468,
                                                                       2504, 9794, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25076, 0, 3,
                                                                       21788, 8258, 21956, 2504,
                                                                       2540, 9902, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25292, 0, 3,
                                                                       21956, 8342, 22124, 2540,
                                                                       2576, 10010, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25508, 0, 3,
                                                                       22124, 8426, 22292, 2576,
                                                                       2612, 10118, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25724, 0, 3,
                                                                       22292, 8510, 22460, 2612,
                                                                       2648, 10226, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25940, 0, 3,
                                                                       22460, 8594, 22628, 2648,
                                                                       2684, 10334, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26156, 0, 3,
                                                                       22628, 8678, 22796, 2684,
                                                                       2720, 10442, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26372, 0, 3,
                                                                       22796, 8762, 22964, 2720,
                                                                       2756, 10550, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26588, 0, 3,
                                                                       23132, 8930, 23348, 2828,
                                                                       2873, 10658, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26858, 0, 3,
                                                                       23348, 9038, 23564, 2873,
                                                                       2918, 10793, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27128, 0, 3,
                                                                       23564, 9146, 23780, 2918,
                                                                       2963, 10928, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27398, 0, 3,
                                                                       23780, 9254, 23996, 2963,
                                                                       3008, 11063, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27668, 0, 3,
                                                                       23996, 9362, 24212, 3008,
                                                                       3053, 11198, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27938, 0, 3,
                                                                       24212, 9470, 24428, 3053,
                                                                       3098, 11333, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28208, 0, 3,
                                                                       24428, 9578, 24644, 3098,
                                                                       3143, 11468, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28478, 0, 3,
                                                                       24860, 9794, 25076, 3233,
                                                                       3278, 11603, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28748, 0, 3,
                                                                       25076, 9902, 25292, 3278,
                                                                       3323, 11738, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29018, 0, 3,
                                                                       25292, 10010, 25508, 3323,
                                                                       3368, 11873, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29288, 0, 3,
                                                                       25508, 10118, 25724, 3368,
                                                                       3413, 12008, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29558, 0, 3,
                                                                       25724, 10226, 25940, 3413,
                                                                       3458, 12143, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29828, 0, 3,
                                                                       25940, 10334, 26156, 3458,
                                                                       3503, 12278, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30098, 0, 3,
                                                                       26156, 10442, 26372, 3503,
                                                                       3548, 12413, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30368, 3, 3638,
                                                                       3641, 12560, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30378, 3, 3641,
                                                                       3644, 12566, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30388, 3, 3644,
                                                                       3647, 12572, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30398, 3, 3647,
                                                                       3650, 12578, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30408, 3, 3650,
                                                                       3653, 12584, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30418, 3, 3653,
                                                                       3656, 12590, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30428, 3, 3656,
                                                                       3659, 12596, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30438, 3, 3659,
                                                                       3662, 12602, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30448, 3, 3662,
                                                                       3665, 12608, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30458, 3, 3665,
                                                                       3668, 12614, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30468, 3, 3668,
                                                                       3671, 12620, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30478, 3, 3671,
                                                                       3674, 12626, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30488, 3, 3674,
                                                                       3677, 12632, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30498, 3, 3683,
                                                                       3686, 12650, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30508, 3, 3686,
                                                                       3689, 12656, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30518, 3, 3689,
                                                                       3692, 12662, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30528, 3, 3692,
                                                                       3695, 12668, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30538, 3, 3695,
                                                                       3698, 12674, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30548, 3, 3698,
                                                                       3701, 12680, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30558, 3, 3701,
                                                                       3704, 12686, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30568, 3, 3704,
                                                                       3707, 12692, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30578, 3, 3707,
                                                                       3710, 12698, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30588, 3, 3710,
                                                                       3713, 12704, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30598, 3, 3713,
                                                                       3716, 12710, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30608, 3, 3716,
                                                                       3719, 12716, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30618, 3, 3719,
                                                                       3722, 12722, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30628, 0, 3,
                                                                       30368, 12560, 30378, 3728,
                                                                       3737, 12764, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30658, 0, 3,
                                                                       30378, 12566, 30388, 3737,
                                                                       3746, 12782, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30688, 0, 3,
                                                                       30388, 12572, 30398, 3746,
                                                                       3755, 12800, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30718, 0, 3,
                                                                       30398, 12578, 30408, 3755,
                                                                       3764, 12818, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30748, 0, 3,
                                                                       30408, 12584, 30418, 3764,
                                                                       3773, 12836, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30778, 0, 3,
                                                                       30418, 12590, 30428, 3773,
                                                                       3782, 12854, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30808, 0, 3,
                                                                       30428, 12596, 30438, 3782,
                                                                       3791, 12872, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30838, 0, 3,
                                                                       30438, 12602, 30448, 3791,
                                                                       3800, 12890, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30868, 0, 3,
                                                                       30448, 12608, 30458, 3800,
                                                                       3809, 12908, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30898, 0, 3,
                                                                       30458, 12614, 30468, 3809,
                                                                       3818, 12926, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30928, 0, 3,
                                                                       30468, 12620, 30478, 3818,
                                                                       3827, 12944, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30958, 0, 3,
                                                                       30478, 12626, 30488, 3827,
                                                                       3836, 12962, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 30988, 0, 3,
                                                                       30498, 12650, 30508, 3854,
                                                                       3863, 13016, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31018, 0, 3,
                                                                       30508, 12656, 30518, 3863,
                                                                       3872, 13034, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31048, 0, 3,
                                                                       30518, 12662, 30528, 3872,
                                                                       3881, 13052, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31078, 0, 3,
                                                                       30528, 12668, 30538, 3881,
                                                                       3890, 13070, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31108, 0, 3,
                                                                       30538, 12674, 30548, 3890,
                                                                       3899, 13088, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31138, 0, 3,
                                                                       30548, 12680, 30558, 3899,
                                                                       3908, 13106, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31168, 0, 3,
                                                                       30558, 12686, 30568, 3908,
                                                                       3917, 13124, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31198, 0, 3,
                                                                       30568, 12692, 30578, 3917,
                                                                       3926, 13142, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31228, 0, 3,
                                                                       30578, 12698, 30588, 3926,
                                                                       3935, 13160, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31258, 0, 3,
                                                                       30588, 12704, 30598, 3935,
                                                                       3944, 13178, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31288, 0, 3,
                                                                       30598, 12710, 30608, 3944,
                                                                       3953, 13196, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31318, 0, 3,
                                                                       30608, 12716, 30618, 3953,
                                                                       3962, 13214, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31348, 0, 3,
                                                                       30628, 12764, 30658, 3980,
                                                                       3998, 13304, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31408, 0, 3,
                                                                       30658, 12782, 30688, 3998,
                                                                       4016, 13340, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31468, 0, 3,
                                                                       30688, 12800, 30718, 4016,
                                                                       4034, 13376, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31528, 0, 3,
                                                                       30718, 12818, 30748, 4034,
                                                                       4052, 13412, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31588, 0, 3,
                                                                       30748, 12836, 30778, 4052,
                                                                       4070, 13448, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31648, 0, 3,
                                                                       30778, 12854, 30808, 4070,
                                                                       4088, 13484, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31708, 0, 3,
                                                                       30808, 12872, 30838, 4088,
                                                                       4106, 13520, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31768, 0, 3,
                                                                       30838, 12890, 30868, 4106,
                                                                       4124, 13556, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31828, 0, 3,
                                                                       30868, 12908, 30898, 4124,
                                                                       4142, 13592, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31888, 0, 3,
                                                                       30898, 12926, 30928, 4142,
                                                                       4160, 13628, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31948, 0, 3,
                                                                       30928, 12944, 30958, 4160,
                                                                       4178, 13664, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32008, 0, 3,
                                                                       30988, 13016, 31018, 4214,
                                                                       4232, 13772, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32068, 0, 3,
                                                                       31018, 13034, 31048, 4232,
                                                                       4250, 13808, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32128, 0, 3,
                                                                       31048, 13052, 31078, 4250,
                                                                       4268, 13844, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32188, 0, 3,
                                                                       31078, 13070, 31108, 4268,
                                                                       4286, 13880, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32248, 0, 3,
                                                                       31108, 13088, 31138, 4286,
                                                                       4304, 13916, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32308, 0, 3,
                                                                       31138, 13106, 31168, 4304,
                                                                       4322, 13952, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32368, 0, 3,
                                                                       31168, 13124, 31198, 4322,
                                                                       4340, 13988, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32428, 0, 3,
                                                                       31198, 13142, 31228, 4340,
                                                                       4358, 14024, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32488, 0, 3,
                                                                       31228, 13160, 31258, 4358,
                                                                       4376, 14060, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32548, 0, 3,
                                                                       31258, 13178, 31288, 4376,
                                                                       4394, 14096, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32608, 0, 3,
                                                                       31288, 13196, 31318, 4394,
                                                                       4412, 14132, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32668, 0, 3,
                                                                       31348, 13304, 31408, 4448,
                                                                       4478, 14288, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32768, 0, 3,
                                                                       31408, 13340, 31468, 4478,
                                                                       4508, 14348, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32868, 0, 3,
                                                                       31468, 13376, 31528, 4508,
                                                                       4538, 14408, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32968, 0, 3,
                                                                       31528, 13412, 31588, 4538,
                                                                       4568, 14468, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33068, 0, 3,
                                                                       31588, 13448, 31648, 4568,
                                                                       4598, 14528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33168, 0, 3,
                                                                       31648, 13484, 31708, 4598,
                                                                       4628, 14588, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33268, 0, 3,
                                                                       31708, 13520, 31768, 4628,
                                                                       4658, 14648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33368, 0, 3,
                                                                       31768, 13556, 31828, 4658,
                                                                       4688, 14708, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33468, 0, 3,
                                                                       31828, 13592, 31888, 4688,
                                                                       4718, 14768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33568, 0, 3,
                                                                       31888, 13628, 31948, 4718,
                                                                       4748, 14828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33668, 0, 3,
                                                                       32008, 13772, 32068, 4808,
                                                                       4838, 15008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33768, 0, 3,
                                                                       32068, 13808, 32128, 4838,
                                                                       4868, 15068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33868, 0, 3,
                                                                       32128, 13844, 32188, 4868,
                                                                       4898, 15128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33968, 0, 3,
                                                                       32188, 13880, 32248, 4898,
                                                                       4928, 15188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 34068, 0, 3,
                                                                       32248, 13916, 32308, 4928,
                                                                       4958, 15248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 34168, 0, 3,
                                                                       32308, 13952, 32368, 4958,
                                                                       4988, 15308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 34268, 0, 3,
                                                                       32368, 13988, 32428, 4988,
                                                                       5018, 15368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 34368, 0, 3,
                                                                       32428, 14024, 32488, 5018,
                                                                       5048, 15428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 34468, 0, 3,
                                                                       32488, 14060, 32548, 5048,
                                                                       5078, 15488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 34568, 0, 3,
                                                                       32548, 14096, 32608, 5078,
                                                                       5108, 15548, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34668, 0, 3,
                                                                       32668, 14288, 32768, 5168,
                                                                       5213, 15788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34818, 0, 3,
                                                                       32768, 14348, 32868, 5213,
                                                                       5258, 15878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34968, 0, 3,
                                                                       32868, 14408, 32968, 5258,
                                                                       5303, 15968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35118, 0, 3,
                                                                       32968, 14468, 33068, 5303,
                                                                       5348, 16058, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35268, 0, 3,
                                                                       33068, 14528, 33168, 5348,
                                                                       5393, 16148, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35418, 0, 3,
                                                                       33168, 14588, 33268, 5393,
                                                                       5438, 16238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35568, 0, 3,
                                                                       33268, 14648, 33368, 5438,
                                                                       5483, 16328, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35718, 0, 3,
                                                                       33368, 14708, 33468, 5483,
                                                                       5528, 16418, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 35868, 0, 3,
                                                                       33468, 14768, 33568, 5528,
                                                                       5573, 16508, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36018, 0, 3,
                                                                       33668, 15008, 33768, 5663,
                                                                       5708, 16778, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36168, 0, 3,
                                                                       33768, 15068, 33868, 5708,
                                                                       5753, 16868, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36318, 0, 3,
                                                                       33868, 15128, 33968, 5753,
                                                                       5798, 16958, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36468, 0, 3,
                                                                       33968, 15188, 34068, 5798,
                                                                       5843, 17048, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36618, 0, 3,
                                                                       34068, 15248, 34168, 5843,
                                                                       5888, 17138, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36768, 0, 3,
                                                                       34168, 15308, 34268, 5888,
                                                                       5933, 17228, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36918, 0, 3,
                                                                       34268, 15368, 34368, 5933,
                                                                       5978, 17318, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37068, 0, 3,
                                                                       34368, 15428, 34468, 5978,
                                                                       6023, 17408, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37218, 0, 3,
                                                                       34468, 15488, 34568, 6023,
                                                                       6068, 17498, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37368, 0, 3,
                                                                       34668, 15788, 34818, 6158,
                                                                       6221, 17840, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37578, 0, 3,
                                                                       34818, 15878, 34968, 6221,
                                                                       6284, 17966, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37788, 0, 3,
                                                                       34968, 15968, 35118, 6284,
                                                                       6347, 18092, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37998, 0, 3,
                                                                       35118, 16058, 35268, 6347,
                                                                       6410, 18218, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 38208, 0, 3,
                                                                       35268, 16148, 35418, 6410,
                                                                       6473, 18344, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 38418, 0, 3,
                                                                       35418, 16238, 35568, 6473,
                                                                       6536, 18470, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 38628, 0, 3,
                                                                       35568, 16328, 35718, 6536,
                                                                       6599, 18596, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 38838, 0, 3,
                                                                       35718, 16418, 35868, 6599,
                                                                       6662, 18722, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39048, 0, 3,
                                                                       36018, 16778, 36168, 6788,
                                                                       6851, 19100, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39258, 0, 3,
                                                                       36168, 16868, 36318, 6851,
                                                                       6914, 19226, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39468, 0, 3,
                                                                       36318, 16958, 36468, 6914,
                                                                       6977, 19352, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39678, 0, 3,
                                                                       36468, 17048, 36618, 6977,
                                                                       7040, 19478, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39888, 0, 3,
                                                                       36618, 17138, 36768, 7040,
                                                                       7103, 19604, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 40098, 0, 3,
                                                                       36768, 17228, 36918, 7103,
                                                                       7166, 19730, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 40308, 0, 3,
                                                                       36918, 17318, 37068, 7166,
                                                                       7229, 19856, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 40518, 0, 3,
                                                                       37068, 17408, 37218, 7229,
                                                                       7292, 19982, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 40728, 0, 3,
                                                                       37368, 17840, 37578, 7418,
                                                                       7502, 20444, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 41008, 0, 3,
                                                                       37578, 17966, 37788, 7502,
                                                                       7586, 20612, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 41288, 0, 3,
                                                                       37788, 18092, 37998, 7586,
                                                                       7670, 20780, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 41568, 0, 3,
                                                                       37998, 18218, 38208, 7670,
                                                                       7754, 20948, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 41848, 0, 3,
                                                                       38208, 18344, 38418, 7754,
                                                                       7838, 21116, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 42128, 0, 3,
                                                                       38418, 18470, 38628, 7838,
                                                                       7922, 21284, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 42408, 0, 3,
                                                                       38628, 18596, 38838, 7922,
                                                                       8006, 21452, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 42688, 0, 3,
                                                                       39048, 19100, 39258, 8174,
                                                                       8258, 21956, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 42968, 0, 3,
                                                                       39258, 19226, 39468, 8258,
                                                                       8342, 22124, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 43248, 0, 3,
                                                                       39468, 19352, 39678, 8342,
                                                                       8426, 22292, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 43528, 0, 3,
                                                                       39678, 19478, 39888, 8426,
                                                                       8510, 22460, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 43808, 0, 3,
                                                                       39888, 19604, 40098, 8510,
                                                                       8594, 22628, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44088, 0, 3,
                                                                       40098, 19730, 40308, 8594,
                                                                       8678, 22796, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44368, 0, 3,
                                                                       40308, 19856, 40518, 8678,
                                                                       8762, 22964, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 44648, 0, 3,
                                                                       40728, 20444, 41008, 8930,
                                                                       9038, 23564, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 45008, 0, 3,
                                                                       41008, 20612, 41288, 9038,
                                                                       9146, 23780, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 45368, 0, 3,
                                                                       41288, 20780, 41568, 9146,
                                                                       9254, 23996, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 45728, 0, 3,
                                                                       41568, 20948, 41848, 9254,
                                                                       9362, 24212, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46088, 0, 3,
                                                                       41848, 21116, 42128, 9362,
                                                                       9470, 24428, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46448, 0, 3,
                                                                       42128, 21284, 42408, 9470,
                                                                       9578, 24644, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46808, 0, 3,
                                                                       42688, 21956, 42968, 9794,
                                                                       9902, 25292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47168, 0, 3,
                                                                       42968, 22124, 43248, 9902,
                                                                       10010, 25508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47528, 0, 3,
                                                                       43248, 22292, 43528,
                                                                       10010, 10118, 25724,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47888, 0, 3,
                                                                       43528, 22460, 43808,
                                                                       10118, 10226, 25940,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 48248, 0, 3,
                                                                       43808, 22628, 44088,
                                                                       10226, 10334, 26156,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 48608, 0, 3,
                                                                       44088, 22796, 44368,
                                                                       10334, 10442, 26372,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 48968, 0, 3,
                                                                       44648, 23564, 45008,
                                                                       10658, 10793, 27128,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 49418, 0, 3,
                                                                       45008, 23780, 45368,
                                                                       10793, 10928, 27398,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 49868, 0, 3,
                                                                       45368, 23996, 45728,
                                                                       10928, 11063, 27668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 50318, 0, 3,
                                                                       45728, 24212, 46088,
                                                                       11063, 11198, 27938,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 50768, 0, 3,
                                                                       46088, 24428, 46448,
                                                                       11198, 11333, 28208,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 51218, 0, 3,
                                                                       46808, 25292, 47168,
                                                                       11603, 11738, 29018,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 51668, 0, 3,
                                                                       47168, 25508, 47528,
                                                                       11738, 11873, 29288,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 52118, 0, 3,
                                                                       47528, 25724, 47888,
                                                                       11873, 12008, 29558,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 52568, 0, 3,
                                                                       47888, 25940, 48248,
                                                                       12008, 12143, 29828,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 53018, 0, 3,
                                                                       48248, 26156, 48608,
                                                                       12143, 12278, 30098,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53468, 3, 12548,
                                                                       12554, 30368, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53483, 3, 12554,
                                                                       12560, 30378, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53498, 3, 12560,
                                                                       12566, 30388, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53513, 3, 12566,
                                                                       12572, 30398, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53528, 3, 12572,
                                                                       12578, 30408, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53543, 3, 12578,
                                                                       12584, 30418, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53558, 3, 12584,
                                                                       12590, 30428, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53573, 3, 12590,
                                                                       12596, 30438, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53588, 3, 12596,
                                                                       12602, 30448, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53603, 3, 12602,
                                                                       12608, 30458, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53618, 3, 12608,
                                                                       12614, 30468, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53633, 3, 12614,
                                                                       12620, 30478, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53648, 3, 12620,
                                                                       12626, 30488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53663, 3, 12638,
                                                                       12644, 30498, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53678, 3, 12644,
                                                                       12650, 30508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53693, 3, 12650,
                                                                       12656, 30518, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53708, 3, 12656,
                                                                       12662, 30528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53723, 3, 12662,
                                                                       12668, 30538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53738, 3, 12668,
                                                                       12674, 30548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53753, 3, 12674,
                                                                       12680, 30558, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53768, 3, 12680,
                                                                       12686, 30568, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53783, 3, 12686,
                                                                       12692, 30578, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53798, 3, 12692,
                                                                       12698, 30588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53813, 3, 12698,
                                                                       12704, 30598, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53828, 3, 12704,
                                                                       12710, 30608, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 53843, 3, 12710,
                                                                       12716, 30618, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53858, 0, 3,
                                                                       53468, 30368, 53483,
                                                                       12728, 12746, 30628,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53903, 0, 3,
                                                                       53483, 30378, 53498,
                                                                       12746, 12764, 30658,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53948, 0, 3,
                                                                       53498, 30388, 53513,
                                                                       12764, 12782, 30688,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53993, 0, 3,
                                                                       53513, 30398, 53528,
                                                                       12782, 12800, 30718,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54038, 0, 3,
                                                                       53528, 30408, 53543,
                                                                       12800, 12818, 30748,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54083, 0, 3,
                                                                       53543, 30418, 53558,
                                                                       12818, 12836, 30778,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54128, 0, 3,
                                                                       53558, 30428, 53573,
                                                                       12836, 12854, 30808,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54173, 0, 3,
                                                                       53573, 30438, 53588,
                                                                       12854, 12872, 30838,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54218, 0, 3,
                                                                       53588, 30448, 53603,
                                                                       12872, 12890, 30868,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54263, 0, 3,
                                                                       53603, 30458, 53618,
                                                                       12890, 12908, 30898,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54308, 0, 3,
                                                                       53618, 30468, 53633,
                                                                       12908, 12926, 30928,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54353, 0, 3,
                                                                       53633, 30478, 53648,
                                                                       12926, 12944, 30958,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54398, 0, 3,
                                                                       53663, 30498, 53678,
                                                                       12980, 12998, 30988,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54443, 0, 3,
                                                                       53678, 30508, 53693,
                                                                       12998, 13016, 31018,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54488, 0, 3,
                                                                       53693, 30518, 53708,
                                                                       13016, 13034, 31048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54533, 0, 3,
                                                                       53708, 30528, 53723,
                                                                       13034, 13052, 31078,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54578, 0, 3,
                                                                       53723, 30538, 53738,
                                                                       13052, 13070, 31108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54623, 0, 3,
                                                                       53738, 30548, 53753,
                                                                       13070, 13088, 31138,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54668, 0, 3,
                                                                       53753, 30558, 53768,
                                                                       13088, 13106, 31168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54713, 0, 3,
                                                                       53768, 30568, 53783,
                                                                       13106, 13124, 31198,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54758, 0, 3,
                                                                       53783, 30578, 53798,
                                                                       13124, 13142, 31228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54803, 0, 3,
                                                                       53798, 30588, 53813,
                                                                       13142, 13160, 31258,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54848, 0, 3,
                                                                       53813, 30598, 53828,
                                                                       13160, 13178, 31288,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 54893, 0, 3,
                                                                       53828, 30608, 53843,
                                                                       13178, 13196, 31318,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54938, 0, 3,
                                                                       53858, 30628, 53903,
                                                                       13232, 13268, 31348,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55028, 0, 3,
                                                                       53903, 30658, 53948,
                                                                       13268, 13304, 31408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55118, 0, 3,
                                                                       53948, 30688, 53993,
                                                                       13304, 13340, 31468,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55208, 0, 3,
                                                                       53993, 30718, 54038,
                                                                       13340, 13376, 31528,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55298, 0, 3,
                                                                       54038, 30748, 54083,
                                                                       13376, 13412, 31588,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55388, 0, 3,
                                                                       54083, 30778, 54128,
                                                                       13412, 13448, 31648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55478, 0, 3,
                                                                       54128, 30808, 54173,
                                                                       13448, 13484, 31708,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55568, 0, 3,
                                                                       54173, 30838, 54218,
                                                                       13484, 13520, 31768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55658, 0, 3,
                                                                       54218, 30868, 54263,
                                                                       13520, 13556, 31828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55748, 0, 3,
                                                                       54263, 30898, 54308,
                                                                       13556, 13592, 31888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55838, 0, 3,
                                                                       54308, 30928, 54353,
                                                                       13592, 13628, 31948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55928, 0, 3,
                                                                       54398, 30988, 54443,
                                                                       13700, 13736, 32008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 56018, 0, 3,
                                                                       54443, 31018, 54488,
                                                                       13736, 13772, 32068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 56108, 0, 3,
                                                                       54488, 31048, 54533,
                                                                       13772, 13808, 32128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 56198, 0, 3,
                                                                       54533, 31078, 54578,
                                                                       13808, 13844, 32188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 56288, 0, 3,
                                                                       54578, 31108, 54623,
                                                                       13844, 13880, 32248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 56378, 0, 3,
                                                                       54623, 31138, 54668,
                                                                       13880, 13916, 32308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 56468, 0, 3,
                                                                       54668, 31168, 54713,
                                                                       13916, 13952, 32368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 56558, 0, 3,
                                                                       54713, 31198, 54758,
                                                                       13952, 13988, 32428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 56648, 0, 3,
                                                                       54758, 31228, 54803,
                                                                       13988, 14024, 32488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 56738, 0, 3,
                                                                       54803, 31258, 54848,
                                                                       14024, 14060, 32548,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 56828, 0, 3,
                                                                       54848, 31288, 54893,
                                                                       14060, 14096, 32608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 56918, 0, 3,
                                                                       54938, 31348, 55028,
                                                                       14168, 14228, 32668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57068, 0, 3,
                                                                       55028, 31408, 55118,
                                                                       14228, 14288, 32768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57218, 0, 3,
                                                                       55118, 31468, 55208,
                                                                       14288, 14348, 32868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57368, 0, 3,
                                                                       55208, 31528, 55298,
                                                                       14348, 14408, 32968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57518, 0, 3,
                                                                       55298, 31588, 55388,
                                                                       14408, 14468, 33068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57668, 0, 3,
                                                                       55388, 31648, 55478,
                                                                       14468, 14528, 33168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57818, 0, 3,
                                                                       55478, 31708, 55568,
                                                                       14528, 14588, 33268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57968, 0, 3,
                                                                       55568, 31768, 55658,
                                                                       14588, 14648, 33368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 58118, 0, 3,
                                                                       55658, 31828, 55748,
                                                                       14648, 14708, 33468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 58268, 0, 3,
                                                                       55748, 31888, 55838,
                                                                       14708, 14768, 33568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 58418, 0, 3,
                                                                       55928, 32008, 56018,
                                                                       14888, 14948, 33668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 58568, 0, 3,
                                                                       56018, 32068, 56108,
                                                                       14948, 15008, 33768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 58718, 0, 3,
                                                                       56108, 32128, 56198,
                                                                       15008, 15068, 33868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 58868, 0, 3,
                                                                       56198, 32188, 56288,
                                                                       15068, 15128, 33968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 59018, 0, 3,
                                                                       56288, 32248, 56378,
                                                                       15128, 15188, 34068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 59168, 0, 3,
                                                                       56378, 32308, 56468,
                                                                       15188, 15248, 34168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 59318, 0, 3,
                                                                       56468, 32368, 56558,
                                                                       15248, 15308, 34268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 59468, 0, 3,
                                                                       56558, 32428, 56648,
                                                                       15308, 15368, 34368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 59618, 0, 3,
                                                                       56648, 32488, 56738,
                                                                       15368, 15428, 34468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 59768, 0, 3,
                                                                       56738, 32548, 56828,
                                                                       15428, 15488, 34568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 59918, 0, 3,
                                                                       56918, 32668, 57068,
                                                                       15608, 15698, 34668,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 60143, 0, 3,
                                                                       57068, 32768, 57218,
                                                                       15698, 15788, 34818,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 60368, 0, 3,
                                                                       57218, 32868, 57368,
                                                                       15788, 15878, 34968,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 60593, 0, 3,
                                                                       57368, 32968, 57518,
                                                                       15878, 15968, 35118,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 60818, 0, 3,
                                                                       57518, 33068, 57668,
                                                                       15968, 16058, 35268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 61043, 0, 3,
                                                                       57668, 33168, 57818,
                                                                       16058, 16148, 35418,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 61268, 0, 3,
                                                                       57818, 33268, 57968,
                                                                       16148, 16238, 35568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 61493, 0, 3,
                                                                       57968, 33368, 58118,
                                                                       16238, 16328, 35718,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 61718, 0, 3,
                                                                       58118, 33468, 58268,
                                                                       16328, 16418, 35868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 61943, 0, 3,
                                                                       58418, 33668, 58568,
                                                                       16598, 16688, 36018,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62168, 0, 3,
                                                                       58568, 33768, 58718,
                                                                       16688, 16778, 36168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62393, 0, 3,
                                                                       58718, 33868, 58868,
                                                                       16778, 16868, 36318,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62618, 0, 3,
                                                                       58868, 33968, 59018,
                                                                       16868, 16958, 36468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62843, 0, 3,
                                                                       59018, 34068, 59168,
                                                                       16958, 17048, 36618,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63068, 0, 3,
                                                                       59168, 34168, 59318,
                                                                       17048, 17138, 36768,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63293, 0, 3,
                                                                       59318, 34268, 59468,
                                                                       17138, 17228, 36918,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63518, 0, 3,
                                                                       59468, 34368, 59618,
                                                                       17228, 17318, 37068,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63743, 0, 3,
                                                                       59618, 34468, 59768,
                                                                       17318, 17408, 37218,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 63968, 0, 3,
                                                                       59918, 34668, 60143,
                                                                       17588, 17714, 37368,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 64283, 0, 3,
                                                                       60143, 34818, 60368,
                                                                       17714, 17840, 37578,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 64598, 0, 3,
                                                                       60368, 34968, 60593,
                                                                       17840, 17966, 37788,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 64913, 0, 3,
                                                                       60593, 35118, 60818,
                                                                       17966, 18092, 37998,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65228, 0, 3,
                                                                       60818, 35268, 61043,
                                                                       18092, 18218, 38208,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65543, 0, 3,
                                                                       61043, 35418, 61268,
                                                                       18218, 18344, 38418,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65858, 0, 3,
                                                                       61268, 35568, 61493,
                                                                       18344, 18470, 38628,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 66173, 0, 3,
                                                                       61493, 35718, 61718,
                                                                       18470, 18596, 38838,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 66488, 0, 3,
                                                                       61943, 36018, 62168,
                                                                       18848, 18974, 39048,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 66803, 0, 3,
                                                                       62168, 36168, 62393,
                                                                       18974, 19100, 39258,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 67118, 0, 3,
                                                                       62393, 36318, 62618,
                                                                       19100, 19226, 39468,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 67433, 0, 3,
                                                                       62618, 36468, 62843,
                                                                       19226, 19352, 39678,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 67748, 0, 3,
                                                                       62843, 36618, 63068,
                                                                       19352, 19478, 39888,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 68063, 0, 3,
                                                                       63068, 36768, 63293,
                                                                       19478, 19604, 40098,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 68378, 0, 3,
                                                                       63293, 36918, 63518,
                                                                       19604, 19730, 40308,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 68693, 0, 3,
                                                                       63518, 37068, 63743,
                                                                       19730, 19856, 40518,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 69008, 0, 3,
                                                                       63968, 37368, 64283,
                                                                       20108, 20276, 40728,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 69428, 0, 3,
                                                                       64283, 37578, 64598,
                                                                       20276, 20444, 41008,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 69848, 0, 3,
                                                                       64598, 37788, 64913,
                                                                       20444, 20612, 41288,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 70268, 0, 3,
                                                                       64913, 37998, 65228,
                                                                       20612, 20780, 41568,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 70688, 0, 3,
                                                                       65228, 38208, 65543,
                                                                       20780, 20948, 41848,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 71108, 0, 3,
                                                                       65543, 38418, 65858,
                                                                       20948, 21116, 42128,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 71528, 0, 3,
                                                                       65858, 38628, 66173,
                                                                       21116, 21284, 42408,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 71948, 0, 3,
                                                                       66488, 39048, 66803,
                                                                       21620, 21788, 42688,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 72368, 0, 3,
                                                                       66803, 39258, 67118,
                                                                       21788, 21956, 42968,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 72788, 0, 3,
                                                                       67118, 39468, 67433,
                                                                       21956, 22124, 43248,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 73208, 0, 3,
                                                                       67433, 39678, 67748,
                                                                       22124, 22292, 43528,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 73628, 0, 3,
                                                                       67748, 39888, 68063,
                                                                       22292, 22460, 43808,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 74048, 0, 3,
                                                                       68063, 40098, 68378,
                                                                       22460, 22628, 44088,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 74468, 0, 3,
                                                                       68378, 40308, 68693,
                                                                       22628, 22796, 44368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 74888, 0, 3,
                                                                       69008, 40728, 69428,
                                                                       23132, 23348, 44648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 75428, 0, 3,
                                                                       69428, 41008, 69848,
                                                                       23348, 23564, 45008,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 75968, 0, 3,
                                                                       69848, 41288, 70268,
                                                                       23564, 23780, 45368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 76508, 0, 3,
                                                                       70268, 41568, 70688,
                                                                       23780, 23996, 45728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 77048, 0, 3,
                                                                       70688, 41848, 71108,
                                                                       23996, 24212, 46088,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 77588, 0, 3,
                                                                       71108, 42128, 71528,
                                                                       24212, 24428, 46448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 78128, 0, 3,
                                                                       71948, 42688, 72368,
                                                                       24860, 25076, 46808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 78668, 0, 3,
                                                                       72368, 42968, 72788,
                                                                       25076, 25292, 47168,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 79208, 0, 3,
                                                                       72788, 43248, 73208,
                                                                       25292, 25508, 47528,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 79748, 0, 3,
                                                                       73208, 43528, 73628,
                                                                       25508, 25724, 47888,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 80288, 0, 3,
                                                                       73628, 43808, 74048,
                                                                       25724, 25940, 48248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 80828, 0, 3,
                                                                       74048, 44088, 74468,
                                                                       25940, 26156, 48608,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 81368, 0, 3,
                                                                       74888, 44648, 75428,
                                                                       26588, 26858, 48968,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 82043, 0, 3,
                                                                       75428, 45008, 75968,
                                                                       26858, 27128, 49418,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 82718, 0, 3,
                                                                       75968, 45368, 76508,
                                                                       27128, 27398, 49868,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 83393, 0, 3,
                                                                       76508, 45728, 77048,
                                                                       27398, 27668, 50318,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 84068, 0, 3,
                                                                       77048, 46088, 77588,
                                                                       27668, 27938, 50768,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 84743, 0, 3,
                                                                       78128, 46808, 78668,
                                                                       28478, 28748, 51218,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 85418, 0, 3,
                                                                       78668, 47168, 79208,
                                                                       28748, 29018, 51668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 86093, 0, 3,
                                                                       79208, 47528, 79748,
                                                                       29018, 29288, 52118,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 86768, 0, 3,
                                                                       79748, 47888, 80288,
                                                                       29288, 29558, 52568,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 87443, 0, 3,
                                                                       80288, 48248, 80828,
                                                                       29558, 29828, 53018,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88118, 3, 30368,
                                                                       30378, 53498, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88139, 3, 30378,
                                                                       30388, 53513, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88160, 3, 30388,
                                                                       30398, 53528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88181, 3, 30398,
                                                                       30408, 53543, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88202, 3, 30408,
                                                                       30418, 53558, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88223, 3, 30418,
                                                                       30428, 53573, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88244, 3, 30428,
                                                                       30438, 53588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88265, 3, 30438,
                                                                       30448, 53603, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88286, 3, 30448,
                                                                       30458, 53618, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88307, 3, 30458,
                                                                       30468, 53633, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88328, 3, 30468,
                                                                       30478, 53648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88349, 3, 30498,
                                                                       30508, 53693, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88370, 3, 30508,
                                                                       30518, 53708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88391, 3, 30518,
                                                                       30528, 53723, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88412, 3, 30528,
                                                                       30538, 53738, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88433, 3, 30538,
                                                                       30548, 53753, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88454, 3, 30548,
                                                                       30558, 53768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88475, 3, 30558,
                                                                       30568, 53783, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88496, 3, 30568,
                                                                       30578, 53798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88517, 3, 30578,
                                                                       30588, 53813, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88538, 3, 30588,
                                                                       30598, 53828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 88559, 3, 30598,
                                                                       30608, 53843, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 88580, 0, 3,
                                                                       88118, 53498, 88139,
                                                                       30628, 30658, 53948,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 88643, 0, 3,
                                                                       88139, 53513, 88160,
                                                                       30658, 30688, 53993,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 88706, 0, 3,
                                                                       88160, 53528, 88181,
                                                                       30688, 30718, 54038,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 88769, 0, 3,
                                                                       88181, 53543, 88202,
                                                                       30718, 30748, 54083,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 88832, 0, 3,
                                                                       88202, 53558, 88223,
                                                                       30748, 30778, 54128,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 88895, 0, 3,
                                                                       88223, 53573, 88244,
                                                                       30778, 30808, 54173,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 88958, 0, 3,
                                                                       88244, 53588, 88265,
                                                                       30808, 30838, 54218,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89021, 0, 3,
                                                                       88265, 53603, 88286,
                                                                       30838, 30868, 54263,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89084, 0, 3,
                                                                       88286, 53618, 88307,
                                                                       30868, 30898, 54308,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89147, 0, 3,
                                                                       88307, 53633, 88328,
                                                                       30898, 30928, 54353,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89210, 0, 3,
                                                                       88349, 53693, 88370,
                                                                       30988, 31018, 54488,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89273, 0, 3,
                                                                       88370, 53708, 88391,
                                                                       31018, 31048, 54533,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89336, 0, 3,
                                                                       88391, 53723, 88412,
                                                                       31048, 31078, 54578,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89399, 0, 3,
                                                                       88412, 53738, 88433,
                                                                       31078, 31108, 54623,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89462, 0, 3,
                                                                       88433, 53753, 88454,
                                                                       31108, 31138, 54668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89525, 0, 3,
                                                                       88454, 53768, 88475,
                                                                       31138, 31168, 54713,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89588, 0, 3,
                                                                       88475, 53783, 88496,
                                                                       31168, 31198, 54758,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89651, 0, 3,
                                                                       88496, 53798, 88517,
                                                                       31198, 31228, 54803,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89714, 0, 3,
                                                                       88517, 53813, 88538,
                                                                       31228, 31258, 54848,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 89777, 0, 3,
                                                                       88538, 53828, 88559,
                                                                       31258, 31288, 54893,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 89840, 0, 3,
                                                                       88580, 53948, 88643,
                                                                       31348, 31408, 55118,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 89966, 0, 3,
                                                                       88643, 53993, 88706,
                                                                       31408, 31468, 55208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 90092, 0, 3,
                                                                       88706, 54038, 88769,
                                                                       31468, 31528, 55298,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 90218, 0, 3,
                                                                       88769, 54083, 88832,
                                                                       31528, 31588, 55388,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 90344, 0, 3,
                                                                       88832, 54128, 88895,
                                                                       31588, 31648, 55478,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 90470, 0, 3,
                                                                       88895, 54173, 88958,
                                                                       31648, 31708, 55568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 90596, 0, 3,
                                                                       88958, 54218, 89021,
                                                                       31708, 31768, 55658,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 90722, 0, 3,
                                                                       89021, 54263, 89084,
                                                                       31768, 31828, 55748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 90848, 0, 3,
                                                                       89084, 54308, 89147,
                                                                       31828, 31888, 55838,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 90974, 0, 3,
                                                                       89210, 54488, 89273,
                                                                       32008, 32068, 56108,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91100, 0, 3,
                                                                       89273, 54533, 89336,
                                                                       32068, 32128, 56198,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91226, 0, 3,
                                                                       89336, 54578, 89399,
                                                                       32128, 32188, 56288,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91352, 0, 3,
                                                                       89399, 54623, 89462,
                                                                       32188, 32248, 56378,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91478, 0, 3,
                                                                       89462, 54668, 89525,
                                                                       32248, 32308, 56468,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91604, 0, 3,
                                                                       89525, 54713, 89588,
                                                                       32308, 32368, 56558,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91730, 0, 3,
                                                                       89588, 54758, 89651,
                                                                       32368, 32428, 56648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91856, 0, 3,
                                                                       89651, 54803, 89714,
                                                                       32428, 32488, 56738,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 91982, 0, 3,
                                                                       89714, 54848, 89777,
                                                                       32488, 32548, 56828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 92108, 0, 3,
                                                                       89840, 55118, 89966,
                                                                       32668, 32768, 57218,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 92318, 0, 3,
                                                                       89966, 55208, 90092,
                                                                       32768, 32868, 57368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 92528, 0, 3,
                                                                       90092, 55298, 90218,
                                                                       32868, 32968, 57518,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 92738, 0, 3,
                                                                       90218, 55388, 90344,
                                                                       32968, 33068, 57668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 92948, 0, 3,
                                                                       90344, 55478, 90470,
                                                                       33068, 33168, 57818,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93158, 0, 3,
                                                                       90470, 55568, 90596,
                                                                       33168, 33268, 57968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93368, 0, 3,
                                                                       90596, 55658, 90722,
                                                                       33268, 33368, 58118,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93578, 0, 3,
                                                                       90722, 55748, 90848,
                                                                       33368, 33468, 58268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93788, 0, 3,
                                                                       90974, 56108, 91100,
                                                                       33668, 33768, 58718,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 93998, 0, 3,
                                                                       91100, 56198, 91226,
                                                                       33768, 33868, 58868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94208, 0, 3,
                                                                       91226, 56288, 91352,
                                                                       33868, 33968, 59018,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94418, 0, 3,
                                                                       91352, 56378, 91478,
                                                                       33968, 34068, 59168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94628, 0, 3,
                                                                       91478, 56468, 91604,
                                                                       34068, 34168, 59318,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 94838, 0, 3,
                                                                       91604, 56558, 91730,
                                                                       34168, 34268, 59468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 95048, 0, 3,
                                                                       91730, 56648, 91856,
                                                                       34268, 34368, 59618,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 95258, 0, 3,
                                                                       91856, 56738, 91982,
                                                                       34368, 34468, 59768,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 95468, 0, 3,
                                                                       92108, 57218, 92318,
                                                                       34668, 34818, 60368,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 95783, 0, 3,
                                                                       92318, 57368, 92528,
                                                                       34818, 34968, 60593,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 96098, 0, 3,
                                                                       92528, 57518, 92738,
                                                                       34968, 35118, 60818,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 96413, 0, 3,
                                                                       92738, 57668, 92948,
                                                                       35118, 35268, 61043,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 96728, 0, 3,
                                                                       92948, 57818, 93158,
                                                                       35268, 35418, 61268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 97043, 0, 3,
                                                                       93158, 57968, 93368,
                                                                       35418, 35568, 61493,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 97358, 0, 3,
                                                                       93368, 58118, 93578,
                                                                       35568, 35718, 61718,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 97673, 0, 3,
                                                                       93788, 58718, 93998,
                                                                       36018, 36168, 62393,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 97988, 0, 3,
                                                                       93998, 58868, 94208,
                                                                       36168, 36318, 62618,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 98303, 0, 3,
                                                                       94208, 59018, 94418,
                                                                       36318, 36468, 62843,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 98618, 0, 3,
                                                                       94418, 59168, 94628,
                                                                       36468, 36618, 63068,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 98933, 0, 3,
                                                                       94628, 59318, 94838,
                                                                       36618, 36768, 63293,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 99248, 0, 3,
                                                                       94838, 59468, 95048,
                                                                       36768, 36918, 63518,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 99563, 0, 3,
                                                                       95048, 59618, 95258,
                                                                       36918, 37068, 63743,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 99878, 0, 3,
                                                                       95468, 60368, 95783,
                                                                       37368, 37578, 64598,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 100319, 0, 3,
                                                                       95783, 60593, 96098,
                                                                       37578, 37788, 64913,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 100760, 0, 3,
                                                                       96098, 60818, 96413,
                                                                       37788, 37998, 65228,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 101201, 0, 3,
                                                                       96413, 61043, 96728,
                                                                       37998, 38208, 65543,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 101642, 0, 3,
                                                                       96728, 61268, 97043,
                                                                       38208, 38418, 65858,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 102083, 0, 3,
                                                                       97043, 61493, 97358,
                                                                       38418, 38628, 66173,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 102524, 0, 3,
                                                                       97673, 62393, 97988,
                                                                       39048, 39258, 67118,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 102965, 0, 3,
                                                                       97988, 62618, 98303,
                                                                       39258, 39468, 67433,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 103406, 0, 3,
                                                                       98303, 62843, 98618,
                                                                       39468, 39678, 67748,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 103847, 0, 3,
                                                                       98618, 63068, 98933,
                                                                       39678, 39888, 68063,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 104288, 0, 3,
                                                                       98933, 63293, 99248,
                                                                       39888, 40098, 68378,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 104729, 0, 3,
                                                                       99248, 63518, 99563,
                                                                       40098, 40308, 68693,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 105170, 0, 3,
                                                                       99878, 64598, 100319,
                                                                       40728, 41008, 69848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 105758, 0, 3,
                                                                       100319, 64913, 100760,
                                                                       41008, 41288, 70268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 106346, 0, 3,
                                                                       100760, 65228, 101201,
                                                                       41288, 41568, 70688,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 106934, 0, 3,
                                                                       101201, 65543, 101642,
                                                                       41568, 41848, 71108,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 107522, 0, 3,
                                                                       101642, 65858, 102083,
                                                                       41848, 42128, 71528,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 108110, 0, 3,
                                                                       102524, 67118, 102965,
                                                                       42688, 42968, 72788,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 108698, 0, 3,
                                                                       102965, 67433, 103406,
                                                                       42968, 43248, 73208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 109286, 0, 3,
                                                                       103406, 67748, 103847,
                                                                       43248, 43528, 73628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 109874, 0, 3,
                                                                       103847, 68063, 104288,
                                                                       43528, 43808, 74048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 110462, 0, 3,
                                                                       104288, 68378, 104729,
                                                                       43808, 44088, 74468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 111050, 0, 3,
                                                                       105170, 69848, 105758,
                                                                       44648, 45008, 75968,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 111806, 0, 3,
                                                                       105758, 70268, 106346,
                                                                       45008, 45368, 76508,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 112562, 0, 3,
                                                                       106346, 70688, 106934,
                                                                       45368, 45728, 77048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 113318, 0, 3,
                                                                       106934, 71108, 107522,
                                                                       45728, 46088, 77588,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 114074, 0, 3,
                                                                       108110, 72788, 108698,
                                                                       46808, 47168, 79208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 114830, 0, 3,
                                                                       108698, 73208, 109286,
                                                                       47168, 47528, 79748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 115586, 0, 3,
                                                                       109286, 73628, 109874,
                                                                       47528, 47888, 80288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 116342, 0, 3,
                                                                       109874, 74048, 110462,
                                                                       47888, 48248, 80828,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 117098, 0, 3,
                                                                       111050, 75968, 111806,
                                                                       48968, 49418, 82718,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 118043, 0, 3,
                                                                       111806, 76508, 112562,
                                                                       49418, 49868, 83393,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 118988, 0, 3,
                                                                       112562, 77048, 113318,
                                                                       49868, 50318, 84068,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 119933, 0, 3,
                                                                       114074, 79208, 114830,
                                                                       51218, 51668, 86093,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 120878, 0, 3,
                                                                       114830, 79748, 115586,
                                                                       51668, 52118, 86768,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 121823, 0, 3,
                                                                       115586, 80288, 116342,
                                                                       52118, 52568, 87443,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 122768, 3, 53468,
                                                                       53483, 88118, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 122796, 3, 53483,
                                                                       53498, 88139, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 122824, 3, 53498,
                                                                       53513, 88160, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 122852, 3, 53513,
                                                                       53528, 88181, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 122880, 3, 53528,
                                                                       53543, 88202, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 122908, 3, 53543,
                                                                       53558, 88223, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 122936, 3, 53558,
                                                                       53573, 88244, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 122964, 3, 53573,
                                                                       53588, 88265, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 122992, 3, 53588,
                                                                       53603, 88286, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123020, 3, 53603,
                                                                       53618, 88307, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123048, 3, 53618,
                                                                       53633, 88328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123076, 3, 53663,
                                                                       53678, 88349, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123104, 3, 53678,
                                                                       53693, 88370, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123132, 3, 53693,
                                                                       53708, 88391, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123160, 3, 53708,
                                                                       53723, 88412, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123188, 3, 53723,
                                                                       53738, 88433, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123216, 3, 53738,
                                                                       53753, 88454, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123244, 3, 53753,
                                                                       53768, 88475, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123272, 3, 53768,
                                                                       53783, 88496, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123300, 3, 53783,
                                                                       53798, 88517, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123328, 3, 53798,
                                                                       53813, 88538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 123356, 3, 53813,
                                                                       53828, 88559, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 123384, 0, 3,
                                                                       122768, 88118, 122796,
                                                                       53858, 53903, 88580,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 123468, 0, 3,
                                                                       122796, 88139, 122824,
                                                                       53903, 53948, 88643,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 123552, 0, 3,
                                                                       122824, 88160, 122852,
                                                                       53948, 53993, 88706,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 123636, 0, 3,
                                                                       122852, 88181, 122880,
                                                                       53993, 54038, 88769,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 123720, 0, 3,
                                                                       122880, 88202, 122908,
                                                                       54038, 54083, 88832,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 123804, 0, 3,
                                                                       122908, 88223, 122936,
                                                                       54083, 54128, 88895,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 123888, 0, 3,
                                                                       122936, 88244, 122964,
                                                                       54128, 54173, 88958,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 123972, 0, 3,
                                                                       122964, 88265, 122992,
                                                                       54173, 54218, 89021,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124056, 0, 3,
                                                                       122992, 88286, 123020,
                                                                       54218, 54263, 89084,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124140, 0, 3,
                                                                       123020, 88307, 123048,
                                                                       54263, 54308, 89147,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124224, 0, 3,
                                                                       123076, 88349, 123104,
                                                                       54398, 54443, 89210,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124308, 0, 3,
                                                                       123104, 88370, 123132,
                                                                       54443, 54488, 89273,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124392, 0, 3,
                                                                       123132, 88391, 123160,
                                                                       54488, 54533, 89336,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124476, 0, 3,
                                                                       123160, 88412, 123188,
                                                                       54533, 54578, 89399,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124560, 0, 3,
                                                                       123188, 88433, 123216,
                                                                       54578, 54623, 89462,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124644, 0, 3,
                                                                       123216, 88454, 123244,
                                                                       54623, 54668, 89525,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124728, 0, 3,
                                                                       123244, 88475, 123272,
                                                                       54668, 54713, 89588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124812, 0, 3,
                                                                       123272, 88496, 123300,
                                                                       54713, 54758, 89651,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124896, 0, 3,
                                                                       123300, 88517, 123328,
                                                                       54758, 54803, 89714,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 124980, 0, 3,
                                                                       123328, 88538, 123356,
                                                                       54803, 54848, 89777,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 125064, 0, 3,
                                                                       123384, 88580, 123468,
                                                                       54938, 55028, 89840,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 125232, 0, 3,
                                                                       123468, 88643, 123552,
                                                                       55028, 55118, 89966,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 125400, 0, 3,
                                                                       123552, 88706, 123636,
                                                                       55118, 55208, 90092,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 125568, 0, 3,
                                                                       123636, 88769, 123720,
                                                                       55208, 55298, 90218,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 125736, 0, 3,
                                                                       123720, 88832, 123804,
                                                                       55298, 55388, 90344,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 125904, 0, 3,
                                                                       123804, 88895, 123888,
                                                                       55388, 55478, 90470,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 126072, 0, 3,
                                                                       123888, 88958, 123972,
                                                                       55478, 55568, 90596,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 126240, 0, 3,
                                                                       123972, 89021, 124056,
                                                                       55568, 55658, 90722,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 126408, 0, 3,
                                                                       124056, 89084, 124140,
                                                                       55658, 55748, 90848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 126576, 0, 3,
                                                                       124224, 89210, 124308,
                                                                       55928, 56018, 90974,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 126744, 0, 3,
                                                                       124308, 89273, 124392,
                                                                       56018, 56108, 91100,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 126912, 0, 3,
                                                                       124392, 89336, 124476,
                                                                       56108, 56198, 91226,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 127080, 0, 3,
                                                                       124476, 89399, 124560,
                                                                       56198, 56288, 91352,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 127248, 0, 3,
                                                                       124560, 89462, 124644,
                                                                       56288, 56378, 91478,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 127416, 0, 3,
                                                                       124644, 89525, 124728,
                                                                       56378, 56468, 91604,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 127584, 0, 3,
                                                                       124728, 89588, 124812,
                                                                       56468, 56558, 91730,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 127752, 0, 3,
                                                                       124812, 89651, 124896,
                                                                       56558, 56648, 91856,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 127920, 0, 3,
                                                                       124896, 89714, 124980,
                                                                       56648, 56738, 91982,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 128088, 0, 3,
                                                                       125064, 89840, 125232,
                                                                       56918, 57068, 92108,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 128368, 0, 3,
                                                                       125232, 89966, 125400,
                                                                       57068, 57218, 92318,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 128648, 0, 3,
                                                                       125400, 90092, 125568,
                                                                       57218, 57368, 92528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 128928, 0, 3,
                                                                       125568, 90218, 125736,
                                                                       57368, 57518, 92738,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 129208, 0, 3,
                                                                       125736, 90344, 125904,
                                                                       57518, 57668, 92948,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 129488, 0, 3,
                                                                       125904, 90470, 126072,
                                                                       57668, 57818, 93158,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 129768, 0, 3,
                                                                       126072, 90596, 126240,
                                                                       57818, 57968, 93368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 130048, 0, 3,
                                                                       126240, 90722, 126408,
                                                                       57968, 58118, 93578,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 130328, 0, 3,
                                                                       126576, 90974, 126744,
                                                                       58418, 58568, 93788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 130608, 0, 3,
                                                                       126744, 91100, 126912,
                                                                       58568, 58718, 93998,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 130888, 0, 3,
                                                                       126912, 91226, 127080,
                                                                       58718, 58868, 94208,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 131168, 0, 3,
                                                                       127080, 91352, 127248,
                                                                       58868, 59018, 94418,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 131448, 0, 3,
                                                                       127248, 91478, 127416,
                                                                       59018, 59168, 94628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 131728, 0, 3,
                                                                       127416, 91604, 127584,
                                                                       59168, 59318, 94838,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 132008, 0, 3,
                                                                       127584, 91730, 127752,
                                                                       59318, 59468, 95048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 132288, 0, 3,
                                                                       127752, 91856, 127920,
                                                                       59468, 59618, 95258,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 132568, 0, 3,
                                                                       128088, 92108, 128368,
                                                                       59918, 60143, 95468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 132988, 0, 3,
                                                                       128368, 92318, 128648,
                                                                       60143, 60368, 95783,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 133408, 0, 3,
                                                                       128648, 92528, 128928,
                                                                       60368, 60593, 96098,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 133828, 0, 3,
                                                                       128928, 92738, 129208,
                                                                       60593, 60818, 96413,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 134248, 0, 3,
                                                                       129208, 92948, 129488,
                                                                       60818, 61043, 96728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 134668, 0, 3,
                                                                       129488, 93158, 129768,
                                                                       61043, 61268, 97043,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 135088, 0, 3,
                                                                       129768, 93368, 130048,
                                                                       61268, 61493, 97358,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 135508, 0, 3,
                                                                       130328, 93788, 130608,
                                                                       61943, 62168, 97673,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 135928, 0, 3,
                                                                       130608, 93998, 130888,
                                                                       62168, 62393, 97988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 136348, 0, 3,
                                                                       130888, 94208, 131168,
                                                                       62393, 62618, 98303,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 136768, 0, 3,
                                                                       131168, 94418, 131448,
                                                                       62618, 62843, 98618,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 137188, 0, 3,
                                                                       131448, 94628, 131728,
                                                                       62843, 63068, 98933,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 137608, 0, 3,
                                                                       131728, 94838, 132008,
                                                                       63068, 63293, 99248,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 138028, 0, 3,
                                                                       132008, 95048, 132288,
                                                                       63293, 63518, 99563,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 138448, 0, 3,
                                                                       132568, 95468, 132988,
                                                                       63968, 64283, 99878,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 139036, 0, 3,
                                                                       132988, 95783, 133408,
                                                                       64283, 64598, 100319,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 139624, 0, 3,
                                                                       133408, 96098, 133828,
                                                                       64598, 64913, 100760,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 140212, 0, 3,
                                                                       133828, 96413, 134248,
                                                                       64913, 65228, 101201,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 140800, 0, 3,
                                                                       134248, 96728, 134668,
                                                                       65228, 65543, 101642,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 141388, 0, 3,
                                                                       134668, 97043, 135088,
                                                                       65543, 65858, 102083,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 141976, 0, 3,
                                                                       135508, 97673, 135928,
                                                                       66488, 66803, 102524,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 142564, 0, 3,
                                                                       135928, 97988, 136348,
                                                                       66803, 67118, 102965,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 143152, 0, 3,
                                                                       136348, 98303, 136768,
                                                                       67118, 67433, 103406,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 143740, 0, 3,
                                                                       136768, 98618, 137188,
                                                                       67433, 67748, 103847,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 144328, 0, 3,
                                                                       137188, 98933, 137608,
                                                                       67748, 68063, 104288,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 144916, 0, 3,
                                                                       137608, 99248, 138028,
                                                                       68063, 68378, 104729,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 145504, 0, 3,
                                                                       138448, 99878, 139036,
                                                                       69008, 69428, 105170,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 146288, 0, 3,
                                                                       139036, 100319, 139624,
                                                                       69428, 69848, 105758,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 147072, 0, 3,
                                                                       139624, 100760, 140212,
                                                                       69848, 70268, 106346,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 147856, 0, 3,
                                                                       140212, 101201, 140800,
                                                                       70268, 70688, 106934,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 148640, 0, 3,
                                                                       140800, 101642, 141388,
                                                                       70688, 71108, 107522,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 149424, 0, 3,
                                                                       141976, 102524, 142564,
                                                                       71948, 72368, 108110,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 150208, 0, 3,
                                                                       142564, 102965, 143152,
                                                                       72368, 72788, 108698,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 150992, 0, 3,
                                                                       143152, 103406, 143740,
                                                                       72788, 73208, 109286,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 151776, 0, 3,
                                                                       143740, 103847, 144328,
                                                                       73208, 73628, 109874,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 152560, 0, 3,
                                                                       144328, 104288, 144916,
                                                                       73628, 74048, 110462,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 153344, 0, 3,
                                                                       145504, 105170, 146288,
                                                                       74888, 75428, 111050,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 154352, 0, 3,
                                                                       146288, 105758, 147072,
                                                                       75428, 75968, 111806,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 155360, 0, 3,
                                                                       147072, 106346, 147856,
                                                                       75968, 76508, 112562,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 156368, 0, 3,
                                                                       147856, 106934, 148640,
                                                                       76508, 77048, 113318,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 157376, 0, 3,
                                                                       149424, 108110, 150208,
                                                                       78128, 78668, 114074,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 158384, 0, 3,
                                                                       150208, 108698, 150992,
                                                                       78668, 79208, 114830,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 159392, 0, 3,
                                                                       150992, 109286, 151776,
                                                                       79208, 79748, 115586,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 160400, 0, 3,
                                                                       151776, 109874, 152560,
                                                                       79748, 80288, 116342,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 161408, 0, 3,
                                                                       153344, 111050, 154352,
                                                                       81368, 82043, 117098,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 162668, 0, 3,
                                                                       154352, 111806, 155360,
                                                                       82043, 82718, 118043,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 163928, 0, 3,
                                                                       155360, 112562, 156368,
                                                                       82718, 83393, 118988,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 165188, 0, 3,
                                                                       157376, 114074, 158384,
                                                                       84743, 85418, 119933,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 166448, 0, 3,
                                                                       158384, 114830, 159392,
                                                                       85418, 86093, 120878,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 167708, 0, 3,
                                                                       159392, 115586, 160400,
                                                                       86093, 86768, 121823,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 168968, 3, 88118,
                                                                       88139, 122824, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169004, 3, 88139,
                                                                       88160, 122852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169040, 3, 88160,
                                                                       88181, 122880, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169076, 3, 88181,
                                                                       88202, 122908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169112, 3, 88202,
                                                                       88223, 122936, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169148, 3, 88223,
                                                                       88244, 122964, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169184, 3, 88244,
                                                                       88265, 122992, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169220, 3, 88265,
                                                                       88286, 123020, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169256, 3, 88286,
                                                                       88307, 123048, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169292, 3, 88349,
                                                                       88370, 123132, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169328, 3, 88370,
                                                                       88391, 123160, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169364, 3, 88391,
                                                                       88412, 123188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169400, 3, 88412,
                                                                       88433, 123216, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169436, 3, 88433,
                                                                       88454, 123244, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169472, 3, 88454,
                                                                       88475, 123272, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169508, 3, 88475,
                                                                       88496, 123300, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169544, 3, 88496,
                                                                       88517, 123328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 169580, 3, 88517,
                                                                       88538, 123356, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 169616, 0, 3,
                                                                       168968, 122824, 169004,
                                                                       88580, 88643, 123552,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 169724, 0, 3,
                                                                       169004, 122852, 169040,
                                                                       88643, 88706, 123636,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 169832, 0, 3,
                                                                       169040, 122880, 169076,
                                                                       88706, 88769, 123720,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 169940, 0, 3,
                                                                       169076, 122908, 169112,
                                                                       88769, 88832, 123804,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 170048, 0, 3,
                                                                       169112, 122936, 169148,
                                                                       88832, 88895, 123888,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 170156, 0, 3,
                                                                       169148, 122964, 169184,
                                                                       88895, 88958, 123972,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 170264, 0, 3,
                                                                       169184, 122992, 169220,
                                                                       88958, 89021, 124056,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 170372, 0, 3,
                                                                       169220, 123020, 169256,
                                                                       89021, 89084, 124140,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 170480, 0, 3,
                                                                       169292, 123132, 169328,
                                                                       89210, 89273, 124392,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 170588, 0, 3,
                                                                       169328, 123160, 169364,
                                                                       89273, 89336, 124476,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 170696, 0, 3,
                                                                       169364, 123188, 169400,
                                                                       89336, 89399, 124560,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 170804, 0, 3,
                                                                       169400, 123216, 169436,
                                                                       89399, 89462, 124644,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 170912, 0, 3,
                                                                       169436, 123244, 169472,
                                                                       89462, 89525, 124728,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 171020, 0, 3,
                                                                       169472, 123272, 169508,
                                                                       89525, 89588, 124812,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 171128, 0, 3,
                                                                       169508, 123300, 169544,
                                                                       89588, 89651, 124896,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 171236, 0, 3,
                                                                       169544, 123328, 169580,
                                                                       89651, 89714, 124980,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 171344, 0, 3,
                                                                       169616, 123552, 169724,
                                                                       89840, 89966, 125400,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 171560, 0, 3,
                                                                       169724, 123636, 169832,
                                                                       89966, 90092, 125568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 171776, 0, 3,
                                                                       169832, 123720, 169940,
                                                                       90092, 90218, 125736,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 171992, 0, 3,
                                                                       169940, 123804, 170048,
                                                                       90218, 90344, 125904,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 172208, 0, 3,
                                                                       170048, 123888, 170156,
                                                                       90344, 90470, 126072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 172424, 0, 3,
                                                                       170156, 123972, 170264,
                                                                       90470, 90596, 126240,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 172640, 0, 3,
                                                                       170264, 124056, 170372,
                                                                       90596, 90722, 126408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 172856, 0, 3,
                                                                       170480, 124392, 170588,
                                                                       90974, 91100, 126912,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 173072, 0, 3,
                                                                       170588, 124476, 170696,
                                                                       91100, 91226, 127080,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 173288, 0, 3,
                                                                       170696, 124560, 170804,
                                                                       91226, 91352, 127248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 173504, 0, 3,
                                                                       170804, 124644, 170912,
                                                                       91352, 91478, 127416,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 173720, 0, 3,
                                                                       170912, 124728, 171020,
                                                                       91478, 91604, 127584,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 173936, 0, 3,
                                                                       171020, 124812, 171128,
                                                                       91604, 91730, 127752,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 174152, 0, 3,
                                                                       171128, 124896, 171236,
                                                                       91730, 91856, 127920,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 174368, 0, 3,
                                                                       171344, 125400, 171560,
                                                                       92108, 92318, 128648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 174728, 0, 3,
                                                                       171560, 125568, 171776,
                                                                       92318, 92528, 128928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 175088, 0, 3,
                                                                       171776, 125736, 171992,
                                                                       92528, 92738, 129208,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 175448, 0, 3,
                                                                       171992, 125904, 172208,
                                                                       92738, 92948, 129488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 175808, 0, 3,
                                                                       172208, 126072, 172424,
                                                                       92948, 93158, 129768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 176168, 0, 3,
                                                                       172424, 126240, 172640,
                                                                       93158, 93368, 130048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 176528, 0, 3,
                                                                       172856, 126912, 173072,
                                                                       93788, 93998, 130888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 176888, 0, 3,
                                                                       173072, 127080, 173288,
                                                                       93998, 94208, 131168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 177248, 0, 3,
                                                                       173288, 127248, 173504,
                                                                       94208, 94418, 131448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 177608, 0, 3,
                                                                       173504, 127416, 173720,
                                                                       94418, 94628, 131728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 177968, 0, 3,
                                                                       173720, 127584, 173936,
                                                                       94628, 94838, 132008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 178328, 0, 3,
                                                                       173936, 127752, 174152,
                                                                       94838, 95048, 132288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 178688, 0, 3,
                                                                       174368, 128648, 174728,
                                                                       95468, 95783, 133408,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 179228, 0, 3,
                                                                       174728, 128928, 175088,
                                                                       95783, 96098, 133828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 179768, 0, 3,
                                                                       175088, 129208, 175448,
                                                                       96098, 96413, 134248,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 180308, 0, 3,
                                                                       175448, 129488, 175808,
                                                                       96413, 96728, 134668,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 180848, 0, 3,
                                                                       175808, 129768, 176168,
                                                                       96728, 97043, 135088,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 181388, 0, 3,
                                                                       176528, 130888, 176888,
                                                                       97673, 97988, 136348,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 181928, 0, 3,
                                                                       176888, 131168, 177248,
                                                                       97988, 98303, 136768,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 182468, 0, 3,
                                                                       177248, 131448, 177608,
                                                                       98303, 98618, 137188,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 183008, 0, 3,
                                                                       177608, 131728, 177968,
                                                                       98618, 98933, 137608,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 183548, 0, 3,
                                                                       177968, 132008, 178328,
                                                                       98933, 99248, 138028,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 184088, 0, 3,
                                                                       178688, 133408, 179228,
                                                                       99878, 100319, 139624,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 184844, 0, 3,
                                                                       179228, 133828, 179768,
                                                                       100319, 100760, 140212,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 185600, 0, 3,
                                                                       179768, 134248, 180308,
                                                                       100760, 101201, 140800,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 186356, 0, 3,
                                                                       180308, 134668, 180848,
                                                                       101201, 101642, 141388,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 187112, 0, 3,
                                                                       181388, 136348, 181928,
                                                                       102524, 102965, 143152,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 187868, 0, 3,
                                                                       181928, 136768, 182468,
                                                                       102965, 103406, 143740,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 188624, 0, 3,
                                                                       182468, 137188, 183008,
                                                                       103406, 103847, 144328,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 189380, 0, 3,
                                                                       183008, 137608, 183548,
                                                                       103847, 104288, 144916,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 190136, 0, 3,
                                                                       184088, 139624, 184844,
                                                                       105170, 105758, 147072,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 191144, 0, 3,
                                                                       184844, 140212, 185600,
                                                                       105758, 106346, 147856,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 192152, 0, 3,
                                                                       185600, 140800, 186356,
                                                                       106346, 106934, 148640,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 193160, 0, 3,
                                                                       187112, 143152, 187868,
                                                                       108110, 108698, 150992,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 194168, 0, 3,
                                                                       187868, 143740, 188624,
                                                                       108698, 109286, 151776,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 195176, 0, 3,
                                                                       188624, 144328, 189380,
                                                                       109286, 109874, 152560,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 196184, 0, 3,
                                                                       190136, 147072, 191144,
                                                                       111050, 111806, 155360,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 197480, 0, 3,
                                                                       191144, 147856, 192152,
                                                                       111806, 112562, 156368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 198776, 0, 3,
                                                                       193160, 150992, 194168,
                                                                       114074, 114830, 159392,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 200072, 0, 3,
                                                                       194168, 151776, 195176,
                                                                       114830, 115586, 160400,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 201368, 0, 3,
                                                                       196184, 155360, 197480,
                                                                       117098, 118043, 163928,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 202988, 0, 3,
                                                                       198776, 159392, 200072,
                                                                       119933, 120878, 167708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 204608, 3, 122768,
                                                                       122796, 168968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 204653, 3, 122796,
                                                                       122824, 169004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 204698, 3, 122824,
                                                                       122852, 169040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 204743, 3, 122852,
                                                                       122880, 169076, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 204788, 3, 122880,
                                                                       122908, 169112, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 204833, 3, 122908,
                                                                       122936, 169148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 204878, 3, 122936,
                                                                       122964, 169184, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 204923, 3, 122964,
                                                                       122992, 169220, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 204968, 3, 122992,
                                                                       123020, 169256, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 205013, 3, 123076,
                                                                       123104, 169292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 205058, 3, 123104,
                                                                       123132, 169328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 205103, 3, 123132,
                                                                       123160, 169364, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 205148, 3, 123160,
                                                                       123188, 169400, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 205193, 3, 123188,
                                                                       123216, 169436, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 205238, 3, 123216,
                                                                       123244, 169472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 205283, 3, 123244,
                                                                       123272, 169508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 205328, 3, 123272,
                                                                       123300, 169544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 205373, 3, 123300,
                                                                       123328, 169580, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 205418, 0, 3,
                                                                       204608, 168968, 204653,
                                                                       123384, 123468, 169616,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 205553, 0, 3,
                                                                       204653, 169004, 204698,
                                                                       123468, 123552, 169724,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 205688, 0, 3,
                                                                       204698, 169040, 204743,
                                                                       123552, 123636, 169832,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 205823, 0, 3,
                                                                       204743, 169076, 204788,
                                                                       123636, 123720, 169940,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 205958, 0, 3,
                                                                       204788, 169112, 204833,
                                                                       123720, 123804, 170048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 206093, 0, 3,
                                                                       204833, 169148, 204878,
                                                                       123804, 123888, 170156,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 206228, 0, 3,
                                                                       204878, 169184, 204923,
                                                                       123888, 123972, 170264,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 206363, 0, 3,
                                                                       204923, 169220, 204968,
                                                                       123972, 124056, 170372,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 206498, 0, 3,
                                                                       205013, 169292, 205058,
                                                                       124224, 124308, 170480,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 206633, 0, 3,
                                                                       205058, 169328, 205103,
                                                                       124308, 124392, 170588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 206768, 0, 3,
                                                                       205103, 169364, 205148,
                                                                       124392, 124476, 170696,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 206903, 0, 3,
                                                                       205148, 169400, 205193,
                                                                       124476, 124560, 170804,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 207038, 0, 3,
                                                                       205193, 169436, 205238,
                                                                       124560, 124644, 170912,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 207173, 0, 3,
                                                                       205238, 169472, 205283,
                                                                       124644, 124728, 171020,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 207308, 0, 3,
                                                                       205283, 169508, 205328,
                                                                       124728, 124812, 171128,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 207443, 0, 3,
                                                                       205328, 169544, 205373,
                                                                       124812, 124896, 171236,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 207578, 0, 3,
                                                                       205418, 169616, 205553,
                                                                       125064, 125232, 171344,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 207848, 0, 3,
                                                                       205553, 169724, 205688,
                                                                       125232, 125400, 171560,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 208118, 0, 3,
                                                                       205688, 169832, 205823,
                                                                       125400, 125568, 171776,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 208388, 0, 3,
                                                                       205823, 169940, 205958,
                                                                       125568, 125736, 171992,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 208658, 0, 3,
                                                                       205958, 170048, 206093,
                                                                       125736, 125904, 172208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 208928, 0, 3,
                                                                       206093, 170156, 206228,
                                                                       125904, 126072, 172424,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 209198, 0, 3,
                                                                       206228, 170264, 206363,
                                                                       126072, 126240, 172640,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 209468, 0, 3,
                                                                       206498, 170480, 206633,
                                                                       126576, 126744, 172856,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 209738, 0, 3,
                                                                       206633, 170588, 206768,
                                                                       126744, 126912, 173072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 210008, 0, 3,
                                                                       206768, 170696, 206903,
                                                                       126912, 127080, 173288,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 210278, 0, 3,
                                                                       206903, 170804, 207038,
                                                                       127080, 127248, 173504,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 210548, 0, 3,
                                                                       207038, 170912, 207173,
                                                                       127248, 127416, 173720,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 210818, 0, 3,
                                                                       207173, 171020, 207308,
                                                                       127416, 127584, 173936,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 211088, 0, 3,
                                                                       207308, 171128, 207443,
                                                                       127584, 127752, 174152,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 211358, 0, 3,
                                                                       207578, 171344, 207848,
                                                                       128088, 128368, 174368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 211808, 0, 3,
                                                                       207848, 171560, 208118,
                                                                       128368, 128648, 174728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 212258, 0, 3,
                                                                       208118, 171776, 208388,
                                                                       128648, 128928, 175088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 212708, 0, 3,
                                                                       208388, 171992, 208658,
                                                                       128928, 129208, 175448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 213158, 0, 3,
                                                                       208658, 172208, 208928,
                                                                       129208, 129488, 175808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 213608, 0, 3,
                                                                       208928, 172424, 209198,
                                                                       129488, 129768, 176168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 214058, 0, 3,
                                                                       209468, 172856, 209738,
                                                                       130328, 130608, 176528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 214508, 0, 3,
                                                                       209738, 173072, 210008,
                                                                       130608, 130888, 176888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 214958, 0, 3,
                                                                       210008, 173288, 210278,
                                                                       130888, 131168, 177248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 215408, 0, 3,
                                                                       210278, 173504, 210548,
                                                                       131168, 131448, 177608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 215858, 0, 3,
                                                                       210548, 173720, 210818,
                                                                       131448, 131728, 177968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 216308, 0, 3,
                                                                       210818, 173936, 211088,
                                                                       131728, 132008, 178328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 216758, 0, 3,
                                                                       211358, 174368, 211808,
                                                                       132568, 132988, 178688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 217433, 0, 3,
                                                                       211808, 174728, 212258,
                                                                       132988, 133408, 179228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 218108, 0, 3,
                                                                       212258, 175088, 212708,
                                                                       133408, 133828, 179768,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 218783, 0, 3,
                                                                       212708, 175448, 213158,
                                                                       133828, 134248, 180308,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 219458, 0, 3,
                                                                       213158, 175808, 213608,
                                                                       134248, 134668, 180848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 220133, 0, 3,
                                                                       214058, 176528, 214508,
                                                                       135508, 135928, 181388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 220808, 0, 3,
                                                                       214508, 176888, 214958,
                                                                       135928, 136348, 181928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 221483, 0, 3,
                                                                       214958, 177248, 215408,
                                                                       136348, 136768, 182468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 222158, 0, 3,
                                                                       215408, 177608, 215858,
                                                                       136768, 137188, 183008,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 222833, 0, 3,
                                                                       215858, 177968, 216308,
                                                                       137188, 137608, 183548,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 223508, 0, 3,
                                                                       216758, 178688, 217433,
                                                                       138448, 139036, 184088,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 224453, 0, 3,
                                                                       217433, 179228, 218108,
                                                                       139036, 139624, 184844,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 225398, 0, 3,
                                                                       218108, 179768, 218783,
                                                                       139624, 140212, 185600,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 226343, 0, 3,
                                                                       218783, 180308, 219458,
                                                                       140212, 140800, 186356,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 227288, 0, 3,
                                                                       220133, 181388, 220808,
                                                                       141976, 142564, 187112,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 228233, 0, 3,
                                                                       220808, 181928, 221483,
                                                                       142564, 143152, 187868,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 229178, 0, 3,
                                                                       221483, 182468, 222158,
                                                                       143152, 143740, 188624,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 230123, 0, 3,
                                                                       222158, 183008, 222833,
                                                                       143740, 144328, 189380,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 231068, 0, 3,
                                                                       223508, 184088, 224453,
                                                                       145504, 146288, 190136,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 232328, 0, 3,
                                                                       224453, 184844, 225398,
                                                                       146288, 147072, 191144,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 233588, 0, 3,
                                                                       225398, 185600, 226343,
                                                                       147072, 147856, 192152,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 234848, 0, 3,
                                                                       227288, 187112, 228233,
                                                                       149424, 150208, 193160,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 236108, 0, 3,
                                                                       228233, 187868, 229178,
                                                                       150208, 150992, 194168,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 237368, 0, 3,
                                                                       229178, 188624, 230123,
                                                                       150992, 151776, 195176,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 238628, 0, 3,
                                                                       231068, 190136, 232328,
                                                                       153344, 154352, 196184,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 240248, 0, 3,
                                                                       232328, 191144, 233588,
                                                                       154352, 155360, 197480,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 241868, 0, 3,
                                                                       234848, 193160, 236108,
                                                                       157376, 158384, 198776,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 243488, 0, 3,
                                                                       236108, 194168, 237368,
                                                                       158384, 159392, 200072,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 245108, 0, 3,
                                                                       238628, 196184, 240248,
                                                                       161408, 162668, 201368,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 247133, 0, 3,
                                                                       241868, 198776, 243488,
                                                                       165188, 166448, 202988,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 249158, 223508, 945, ncols);

                    simdfunc::contract_primitives(buffer, 250460, 227288, 945, ncols);

                    simdfunc::contract_primitives(buffer, 251762, 231068, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 253498, 234848, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 255234, 238628, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 257466, 241868, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 259698, 245108, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 262488, 247133, 2025, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 250103, 249158, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 251405, 250460, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 253022, 251762, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 254758, 253498, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 256854, 255234, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 259086, 257466, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 261723, 259698, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 264513, 262488, 45, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 265278, 250103, 253022, 17,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 266349, 251405, 254758, 17,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 267420, 253022, 256854, 17,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 268848, 254758, 259086, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 270276, 256854, 261723, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 272112, 259086, 264513, 17,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 273948, 265278, 267420, 17,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 276090, 266349, 268848, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 278232, 267420, 270276, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 281088, 268848, 272112, 17,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 283944, 273948, 278232, 17,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 287514, 276090, 281088, 17,
                                             nmax);

        simdtrf::transform_f_inner(buffer, 291084, 287514, 21, 17, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 291084, 119, nmax);

        simdtrf::transform_f_inner(buffer, 291084, 283944, 21, 17, nmax);

        simdtrf::transform_h_outer(values + 1309 * nvalues + n * npairs, nvalues, buffer, 291084,
                                   119, nmax);
    }

    for (size_t m = 0; m < 2618; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
