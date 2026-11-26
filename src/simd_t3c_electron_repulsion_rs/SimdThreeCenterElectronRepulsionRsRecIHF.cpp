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


#include "SimdThreeCenterElectronRepulsionRsRecIHF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIH.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKG.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLF.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMD.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransferNP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_ihf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ihf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 133428, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2002 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 133428, 60432, 9926, dimensions);

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

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3198, 0, 3, 2108,
                                                                       2153, 2648, 2703, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3264, 0, 3, 2153,
                                                                       2198, 2703, 2758, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3330, 0, 3, 2198,
                                                                       2243, 2758, 2813, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3396, 0, 3, 2243,
                                                                       2288, 2813, 2868, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3462, 0, 3, 2378,
                                                                       2423, 2923, 2978, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3528, 0, 3, 2423,
                                                                       2468, 2978, 3033, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3594, 0, 3, 2468,
                                                                       2513, 3033, 3088, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3660, 0, 3, 2513,
                                                                       2558, 3088, 3143, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3726, 0, 3, 2648,
                                                                       2703, 3198, 3264, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3804, 0, 3, 2703,
                                                                       2758, 3264, 3330, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3882, 0, 3, 2758,
                                                                       2813, 3330, 3396, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3960, 0, 3, 2923,
                                                                       2978, 3462, 3528, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4038, 0, 3, 2978,
                                                                       3033, 3528, 3594, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4116, 0, 3, 3033,
                                                                       3088, 3594, 3660, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4194, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4197, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4200, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4203, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4206, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4209, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4212, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4215, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4218, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4221, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4224, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4227, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4230, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4233, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4236, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4239, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4242, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4245, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4248, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4251, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4254, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4257, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4260, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4263, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4266, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4269, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4272, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4275, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4278, 3, 7, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4287, 3, 8, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4296, 3, 9, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4305, 3, 10, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4314, 3, 11, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4323, 3, 12, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4332, 3, 13, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4341, 3, 14, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4350, 3, 15, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4359, 3, 16, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4368, 3, 17, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4377, 3, 18, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4386, 3, 19, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4395, 3, 22, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4404, 3, 23, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4413, 3, 24, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4422, 3, 25, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4431, 3, 26, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4440, 3, 27, 90,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4449, 3, 28, 93,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4458, 3, 29, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4467, 3, 30, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4476, 3, 31, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4485, 3, 32, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4494, 3, 33, 108,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4503, 3, 34, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4512, 3, 36, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4530, 3, 39, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4548, 3, 42, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4566, 3, 45, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4584, 3, 48, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4602, 3, 51, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4620, 3, 54, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4638, 3, 57, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4656, 3, 60, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4674, 3, 63, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4692, 3, 66, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4710, 3, 69, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4728, 3, 75, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4746, 3, 78, 192,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4764, 3, 81, 198,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4782, 3, 84, 204,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4800, 3, 87, 210,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4818, 3, 90, 216,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4836, 3, 93, 222,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4854, 3, 96, 228,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4872, 3, 99, 234,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4890, 3, 102, 240,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4908, 3, 105, 246,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4926, 3, 108, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4944, 3, 114, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4974, 3, 120, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5004, 3, 126, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5034, 3, 132, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5064, 3, 138, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5094, 3, 144, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5124, 3, 150, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5154, 3, 156, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5184, 3, 162, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5214, 3, 168, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5244, 3, 174, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5274, 3, 186, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5304, 3, 192, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5334, 3, 198, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5364, 3, 204, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5394, 3, 210, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5424, 3, 216, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5454, 3, 222, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5484, 3, 228, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5514, 3, 234, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5544, 3, 240, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5574, 3, 246, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5604, 3, 258, 478,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5649, 3, 268, 493,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5694, 3, 278, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5739, 3, 288, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5784, 3, 298, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5829, 3, 308, 553,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5874, 3, 318, 568,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5919, 3, 328, 583,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5964, 3, 338, 598,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6009, 3, 348, 613,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6054, 3, 368, 628,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6099, 3, 378, 643,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6144, 3, 388, 658,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6189, 3, 398, 673,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6234, 3, 408, 688,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6279, 3, 418, 703,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6324, 3, 428, 718,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6369, 3, 438, 733,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6414, 3, 448, 748,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6459, 3, 458, 763,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6504, 3, 478, 778,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6567, 3, 493, 799,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6630, 3, 508, 820,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6693, 3, 523, 841,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6756, 3, 538, 862,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6819, 3, 553, 883,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6882, 3, 568, 904,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6945, 3, 583, 925,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7008, 3, 598, 946,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7071, 3, 628, 967,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7134, 3, 643, 988,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7197, 3, 658,
                                                                       1009, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7260, 3, 673,
                                                                       1030, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7323, 3, 688,
                                                                       1051, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7386, 3, 703,
                                                                       1072, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7449, 3, 718,
                                                                       1093, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7512, 3, 733,
                                                                       1114, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7575, 3, 748,
                                                                       1135, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7638, 3, 778,
                                                                       1156, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7722, 3, 799,
                                                                       1184, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7806, 3, 820,
                                                                       1212, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7890, 3, 841,
                                                                       1240, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7974, 3, 862,
                                                                       1268, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8058, 3, 883,
                                                                       1296, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8142, 3, 904,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8226, 3, 925,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8310, 3, 967,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8394, 3, 988,
                                                                       1408, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8478, 3, 1009,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8562, 3, 1030,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8646, 3, 1051,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8730, 3, 1072,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8814, 3, 1093,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8898, 3, 1114,
                                                                       1576, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8982, 3, 1156,
                                                                       1604, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9090, 3, 1184,
                                                                       1640, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9198, 3, 1212,
                                                                       1676, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9306, 3, 1240,
                                                                       1712, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9414, 3, 1268,
                                                                       1748, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9522, 3, 1296,
                                                                       1784, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9630, 3, 1324,
                                                                       1820, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9738, 3, 1380,
                                                                       1856, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9846, 3, 1408,
                                                                       1892, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9954, 3, 1436,
                                                                       1928, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10062, 3, 1464,
                                                                       1964, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10170, 3, 1492,
                                                                       2000, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10278, 3, 1520,
                                                                       2036, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10386, 3, 1548,
                                                                       2072, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10494, 3, 1604,
                                                                       2108, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10629, 3, 1640,
                                                                       2153, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10764, 3, 1676,
                                                                       2198, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10899, 3, 1712,
                                                                       2243, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11034, 3, 1748,
                                                                       2288, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11169, 3, 1784,
                                                                       2333, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11304, 3, 1856,
                                                                       2378, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11439, 3, 1892,
                                                                       2423, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11574, 3, 1928,
                                                                       2468, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11709, 3, 1964,
                                                                       2513, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11844, 3, 2000,
                                                                       2558, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11979, 3, 2036,
                                                                       2603, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12114, 3, 2108,
                                                                       2648, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12279, 3, 2153,
                                                                       2703, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12444, 3, 2198,
                                                                       2758, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12609, 3, 2243,
                                                                       2813, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12774, 3, 2288,
                                                                       2868, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12939, 3, 2378,
                                                                       2923, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13104, 3, 2423,
                                                                       2978, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13269, 3, 2468,
                                                                       3033, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13434, 3, 2513,
                                                                       3088, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13599, 3, 2558,
                                                                       3143, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13764, 3, 2648,
                                                                       3198, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13962, 3, 2703,
                                                                       3264, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14160, 3, 2758,
                                                                       3330, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14358, 3, 2813,
                                                                       3396, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14556, 3, 2923,
                                                                       3462, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14754, 3, 2978,
                                                                       3528, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14952, 3, 3033,
                                                                       3594, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15150, 3, 3088,
                                                                       3660, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15348, 3, 3198,
                                                                       3726, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15582, 3, 3264,
                                                                       3804, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15816, 3, 3330,
                                                                       3882, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16050, 3, 3462,
                                                                       3960, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16284, 3, 3528,
                                                                       4038, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16518, 3, 3594,
                                                                       4116, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16752, 3, 7, 8,
                                                                       4200, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16758, 3, 8, 9,
                                                                       4203, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16764, 3, 9, 10,
                                                                       4206, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16770, 3, 10, 11,
                                                                       4209, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16776, 3, 11, 12,
                                                                       4212, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16782, 3, 12, 13,
                                                                       4215, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16788, 3, 13, 14,
                                                                       4218, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16794, 3, 14, 15,
                                                                       4221, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16800, 3, 15, 16,
                                                                       4224, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16806, 3, 16, 17,
                                                                       4227, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16812, 3, 17, 18,
                                                                       4230, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16818, 3, 18, 19,
                                                                       4233, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16824, 3, 22, 23,
                                                                       4242, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16830, 3, 23, 24,
                                                                       4245, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16836, 3, 24, 25,
                                                                       4248, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16842, 3, 25, 26,
                                                                       4251, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16848, 3, 26, 27,
                                                                       4254, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16854, 3, 27, 28,
                                                                       4257, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16860, 3, 28, 29,
                                                                       4260, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16866, 3, 29, 30,
                                                                       4263, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16872, 3, 30, 31,
                                                                       4266, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16878, 3, 31, 32,
                                                                       4269, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16884, 3, 32, 33,
                                                                       4272, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16890, 3, 33, 34,
                                                                       4275, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16896, 0, 3,
                                                                       16752, 4200, 16758, 4296,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16914, 0, 3,
                                                                       16758, 4203, 16764, 4305,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16932, 0, 3,
                                                                       16764, 4206, 16770, 4314,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16950, 0, 3,
                                                                       16770, 4209, 16776, 4323,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16968, 0, 3,
                                                                       16776, 4212, 16782, 4332,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16986, 0, 3,
                                                                       16782, 4215, 16788, 4341,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17004, 0, 3,
                                                                       16788, 4218, 16794, 4350,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17022, 0, 3,
                                                                       16794, 4221, 16800, 4359,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17040, 0, 3,
                                                                       16800, 4224, 16806, 4368,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17058, 0, 3,
                                                                       16806, 4227, 16812, 4377,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17076, 0, 3,
                                                                       16812, 4230, 16818, 4386,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17094, 0, 3,
                                                                       16824, 4242, 16830, 4413,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17112, 0, 3,
                                                                       16830, 4245, 16836, 4422,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17130, 0, 3,
                                                                       16836, 4248, 16842, 4431,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17148, 0, 3,
                                                                       16842, 4251, 16848, 4440,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17166, 0, 3,
                                                                       16848, 4254, 16854, 4449,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17184, 0, 3,
                                                                       16854, 4257, 16860, 4458,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17202, 0, 3,
                                                                       16860, 4260, 16866, 4467,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17220, 0, 3,
                                                                       16866, 4263, 16872, 4476,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17238, 0, 3,
                                                                       16872, 4266, 16878, 4485,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17256, 0, 3,
                                                                       16878, 4269, 16884, 4494,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17274, 0, 3,
                                                                       16884, 4272, 16890, 4503,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17292, 0, 3,
                                                                       16896, 4296, 16914, 114,
                                                                       120, 4548, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17328, 0, 3,
                                                                       16914, 4305, 16932, 120,
                                                                       126, 4566, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17364, 0, 3,
                                                                       16932, 4314, 16950, 126,
                                                                       132, 4584, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17400, 0, 3,
                                                                       16950, 4323, 16968, 132,
                                                                       138, 4602, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17436, 0, 3,
                                                                       16968, 4332, 16986, 138,
                                                                       144, 4620, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17472, 0, 3,
                                                                       16986, 4341, 17004, 144,
                                                                       150, 4638, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17508, 0, 3,
                                                                       17004, 4350, 17022, 150,
                                                                       156, 4656, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17544, 0, 3,
                                                                       17022, 4359, 17040, 156,
                                                                       162, 4674, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17580, 0, 3,
                                                                       17040, 4368, 17058, 162,
                                                                       168, 4692, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17616, 0, 3,
                                                                       17058, 4377, 17076, 168,
                                                                       174, 4710, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17652, 0, 3,
                                                                       17094, 4413, 17112, 186,
                                                                       192, 4764, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17688, 0, 3,
                                                                       17112, 4422, 17130, 192,
                                                                       198, 4782, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17724, 0, 3,
                                                                       17130, 4431, 17148, 198,
                                                                       204, 4800, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17760, 0, 3,
                                                                       17148, 4440, 17166, 204,
                                                                       210, 4818, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17796, 0, 3,
                                                                       17166, 4449, 17184, 210,
                                                                       216, 4836, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17832, 0, 3,
                                                                       17184, 4458, 17202, 216,
                                                                       222, 4854, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17868, 0, 3,
                                                                       17202, 4467, 17220, 222,
                                                                       228, 4872, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17904, 0, 3,
                                                                       17220, 4476, 17238, 228,
                                                                       234, 4890, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17940, 0, 3,
                                                                       17238, 4485, 17256, 234,
                                                                       240, 4908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17976, 0, 3,
                                                                       17256, 4494, 17274, 240,
                                                                       246, 4926, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18012, 0, 3,
                                                                       17292, 4548, 17328, 258,
                                                                       268, 5004, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18072, 0, 3,
                                                                       17328, 4566, 17364, 268,
                                                                       278, 5034, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18132, 0, 3,
                                                                       17364, 4584, 17400, 278,
                                                                       288, 5064, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18192, 0, 3,
                                                                       17400, 4602, 17436, 288,
                                                                       298, 5094, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18252, 0, 3,
                                                                       17436, 4620, 17472, 298,
                                                                       308, 5124, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18312, 0, 3,
                                                                       17472, 4638, 17508, 308,
                                                                       318, 5154, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18372, 0, 3,
                                                                       17508, 4656, 17544, 318,
                                                                       328, 5184, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18432, 0, 3,
                                                                       17544, 4674, 17580, 328,
                                                                       338, 5214, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18492, 0, 3,
                                                                       17580, 4692, 17616, 338,
                                                                       348, 5244, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18552, 0, 3,
                                                                       17652, 4764, 17688, 368,
                                                                       378, 5334, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18612, 0, 3,
                                                                       17688, 4782, 17724, 378,
                                                                       388, 5364, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18672, 0, 3,
                                                                       17724, 4800, 17760, 388,
                                                                       398, 5394, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18732, 0, 3,
                                                                       17760, 4818, 17796, 398,
                                                                       408, 5424, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18792, 0, 3,
                                                                       17796, 4836, 17832, 408,
                                                                       418, 5454, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18852, 0, 3,
                                                                       17832, 4854, 17868, 418,
                                                                       428, 5484, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18912, 0, 3,
                                                                       17868, 4872, 17904, 428,
                                                                       438, 5514, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18972, 0, 3,
                                                                       17904, 4890, 17940, 438,
                                                                       448, 5544, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19032, 0, 3,
                                                                       17940, 4908, 17976, 448,
                                                                       458, 5574, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19092, 0, 3,
                                                                       18012, 5004, 18072, 478,
                                                                       493, 5694, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19182, 0, 3,
                                                                       18072, 5034, 18132, 493,
                                                                       508, 5739, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19272, 0, 3,
                                                                       18132, 5064, 18192, 508,
                                                                       523, 5784, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19362, 0, 3,
                                                                       18192, 5094, 18252, 523,
                                                                       538, 5829, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19452, 0, 3,
                                                                       18252, 5124, 18312, 538,
                                                                       553, 5874, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19542, 0, 3,
                                                                       18312, 5154, 18372, 553,
                                                                       568, 5919, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19632, 0, 3,
                                                                       18372, 5184, 18432, 568,
                                                                       583, 5964, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19722, 0, 3,
                                                                       18432, 5214, 18492, 583,
                                                                       598, 6009, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19812, 0, 3,
                                                                       18552, 5334, 18612, 628,
                                                                       643, 6144, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19902, 0, 3,
                                                                       18612, 5364, 18672, 643,
                                                                       658, 6189, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19992, 0, 3,
                                                                       18672, 5394, 18732, 658,
                                                                       673, 6234, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20082, 0, 3,
                                                                       18732, 5424, 18792, 673,
                                                                       688, 6279, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20172, 0, 3,
                                                                       18792, 5454, 18852, 688,
                                                                       703, 6324, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20262, 0, 3,
                                                                       18852, 5484, 18912, 703,
                                                                       718, 6369, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20352, 0, 3,
                                                                       18912, 5514, 18972, 718,
                                                                       733, 6414, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20442, 0, 3,
                                                                       18972, 5544, 19032, 733,
                                                                       748, 6459, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 20532, 0, 3,
                                                                       19092, 5694, 19182, 778,
                                                                       799, 6630, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 20658, 0, 3,
                                                                       19182, 5739, 19272, 799,
                                                                       820, 6693, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 20784, 0, 3,
                                                                       19272, 5784, 19362, 820,
                                                                       841, 6756, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 20910, 0, 3,
                                                                       19362, 5829, 19452, 841,
                                                                       862, 6819, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21036, 0, 3,
                                                                       19452, 5874, 19542, 862,
                                                                       883, 6882, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21162, 0, 3,
                                                                       19542, 5919, 19632, 883,
                                                                       904, 6945, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21288, 0, 3,
                                                                       19632, 5964, 19722, 904,
                                                                       925, 7008, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21414, 0, 3,
                                                                       19812, 6144, 19902, 967,
                                                                       988, 7197, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21540, 0, 3,
                                                                       19902, 6189, 19992, 988,
                                                                       1009, 7260, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21666, 0, 3,
                                                                       19992, 6234, 20082, 1009,
                                                                       1030, 7323, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21792, 0, 3,
                                                                       20082, 6279, 20172, 1030,
                                                                       1051, 7386, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21918, 0, 3,
                                                                       20172, 6324, 20262, 1051,
                                                                       1072, 7449, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22044, 0, 3,
                                                                       20262, 6369, 20352, 1072,
                                                                       1093, 7512, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22170, 0, 3,
                                                                       20352, 6414, 20442, 1093,
                                                                       1114, 7575, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22296, 0, 3,
                                                                       20532, 6630, 20658, 1156,
                                                                       1184, 7806, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22464, 0, 3,
                                                                       20658, 6693, 20784, 1184,
                                                                       1212, 7890, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22632, 0, 3,
                                                                       20784, 6756, 20910, 1212,
                                                                       1240, 7974, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22800, 0, 3,
                                                                       20910, 6819, 21036, 1240,
                                                                       1268, 8058, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22968, 0, 3,
                                                                       21036, 6882, 21162, 1268,
                                                                       1296, 8142, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23136, 0, 3,
                                                                       21162, 6945, 21288, 1296,
                                                                       1324, 8226, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23304, 0, 3,
                                                                       21414, 7197, 21540, 1380,
                                                                       1408, 8478, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23472, 0, 3,
                                                                       21540, 7260, 21666, 1408,
                                                                       1436, 8562, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23640, 0, 3,
                                                                       21666, 7323, 21792, 1436,
                                                                       1464, 8646, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23808, 0, 3,
                                                                       21792, 7386, 21918, 1464,
                                                                       1492, 8730, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23976, 0, 3,
                                                                       21918, 7449, 22044, 1492,
                                                                       1520, 8814, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24144, 0, 3,
                                                                       22044, 7512, 22170, 1520,
                                                                       1548, 8898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24312, 0, 3,
                                                                       22296, 7806, 22464, 1604,
                                                                       1640, 9198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24528, 0, 3,
                                                                       22464, 7890, 22632, 1640,
                                                                       1676, 9306, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24744, 0, 3,
                                                                       22632, 7974, 22800, 1676,
                                                                       1712, 9414, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24960, 0, 3,
                                                                       22800, 8058, 22968, 1712,
                                                                       1748, 9522, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25176, 0, 3,
                                                                       22968, 8142, 23136, 1748,
                                                                       1784, 9630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25392, 0, 3,
                                                                       23304, 8478, 23472, 1856,
                                                                       1892, 9954, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25608, 0, 3,
                                                                       23472, 8562, 23640, 1892,
                                                                       1928, 10062, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25824, 0, 3,
                                                                       23640, 8646, 23808, 1928,
                                                                       1964, 10170, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26040, 0, 3,
                                                                       23808, 8730, 23976, 1964,
                                                                       2000, 10278, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26256, 0, 3,
                                                                       23976, 8814, 24144, 2000,
                                                                       2036, 10386, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26472, 0, 3,
                                                                       24312, 9198, 24528, 2108,
                                                                       2153, 10764, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26742, 0, 3,
                                                                       24528, 9306, 24744, 2153,
                                                                       2198, 10899, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27012, 0, 3,
                                                                       24744, 9414, 24960, 2198,
                                                                       2243, 11034, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27282, 0, 3,
                                                                       24960, 9522, 25176, 2243,
                                                                       2288, 11169, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27552, 0, 3,
                                                                       25392, 9954, 25608, 2378,
                                                                       2423, 11574, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27822, 0, 3,
                                                                       25608, 10062, 25824, 2423,
                                                                       2468, 11709, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28092, 0, 3,
                                                                       25824, 10170, 26040, 2468,
                                                                       2513, 11844, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28362, 0, 3,
                                                                       26040, 10278, 26256, 2513,
                                                                       2558, 11979, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 28632, 0, 3,
                                                                       26472, 10764, 26742, 2648,
                                                                       2703, 12444, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 28962, 0, 3,
                                                                       26742, 10899, 27012, 2703,
                                                                       2758, 12609, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29292, 0, 3,
                                                                       27012, 11034, 27282, 2758,
                                                                       2813, 12774, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29622, 0, 3,
                                                                       27552, 11574, 27822, 2923,
                                                                       2978, 13269, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29952, 0, 3,
                                                                       27822, 11709, 28092, 2978,
                                                                       3033, 13434, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 30282, 0, 3,
                                                                       28092, 11844, 28362, 3033,
                                                                       3088, 13599, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 30612, 0, 3,
                                                                       28632, 12444, 28962, 3198,
                                                                       3264, 14160, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 31008, 0, 3,
                                                                       28962, 12609, 29292, 3264,
                                                                       3330, 14358, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 31404, 0, 3,
                                                                       29622, 13269, 29952, 3462,
                                                                       3528, 14952, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 31800, 0, 3,
                                                                       29952, 13434, 30282, 3528,
                                                                       3594, 15150, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 32196, 0, 3,
                                                                       30612, 14160, 31008, 3726,
                                                                       3804, 15816, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 32664, 0, 3,
                                                                       31404, 14952, 31800, 3960,
                                                                       4038, 16518, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33132, 3, 4194,
                                                                       4197, 16752, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33142, 3, 4197,
                                                                       4200, 16758, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33152, 3, 4200,
                                                                       4203, 16764, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33162, 3, 4203,
                                                                       4206, 16770, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33172, 3, 4206,
                                                                       4209, 16776, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33182, 3, 4209,
                                                                       4212, 16782, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33192, 3, 4212,
                                                                       4215, 16788, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33202, 3, 4215,
                                                                       4218, 16794, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33212, 3, 4218,
                                                                       4221, 16800, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33222, 3, 4221,
                                                                       4224, 16806, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33232, 3, 4224,
                                                                       4227, 16812, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33242, 3, 4227,
                                                                       4230, 16818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33252, 3, 4236,
                                                                       4239, 16824, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33262, 3, 4239,
                                                                       4242, 16830, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33272, 3, 4242,
                                                                       4245, 16836, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33282, 3, 4245,
                                                                       4248, 16842, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33292, 3, 4248,
                                                                       4251, 16848, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33302, 3, 4251,
                                                                       4254, 16854, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33312, 3, 4254,
                                                                       4257, 16860, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33322, 3, 4257,
                                                                       4260, 16866, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33332, 3, 4260,
                                                                       4263, 16872, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33342, 3, 4263,
                                                                       4266, 16878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33352, 3, 4266,
                                                                       4269, 16884, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33362, 3, 4269,
                                                                       4272, 16890, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33372, 0, 3,
                                                                       33132, 16752, 33142, 4278,
                                                                       4287, 16896, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33402, 0, 3,
                                                                       33142, 16758, 33152, 4287,
                                                                       4296, 16914, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33432, 0, 3,
                                                                       33152, 16764, 33162, 4296,
                                                                       4305, 16932, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33462, 0, 3,
                                                                       33162, 16770, 33172, 4305,
                                                                       4314, 16950, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33492, 0, 3,
                                                                       33172, 16776, 33182, 4314,
                                                                       4323, 16968, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33522, 0, 3,
                                                                       33182, 16782, 33192, 4323,
                                                                       4332, 16986, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33552, 0, 3,
                                                                       33192, 16788, 33202, 4332,
                                                                       4341, 17004, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33582, 0, 3,
                                                                       33202, 16794, 33212, 4341,
                                                                       4350, 17022, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33612, 0, 3,
                                                                       33212, 16800, 33222, 4350,
                                                                       4359, 17040, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33642, 0, 3,
                                                                       33222, 16806, 33232, 4359,
                                                                       4368, 17058, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33672, 0, 3,
                                                                       33232, 16812, 33242, 4368,
                                                                       4377, 17076, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33702, 0, 3,
                                                                       33252, 16824, 33262, 4395,
                                                                       4404, 17094, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33732, 0, 3,
                                                                       33262, 16830, 33272, 4404,
                                                                       4413, 17112, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33762, 0, 3,
                                                                       33272, 16836, 33282, 4413,
                                                                       4422, 17130, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33792, 0, 3,
                                                                       33282, 16842, 33292, 4422,
                                                                       4431, 17148, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33822, 0, 3,
                                                                       33292, 16848, 33302, 4431,
                                                                       4440, 17166, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33852, 0, 3,
                                                                       33302, 16854, 33312, 4440,
                                                                       4449, 17184, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33882, 0, 3,
                                                                       33312, 16860, 33322, 4449,
                                                                       4458, 17202, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33912, 0, 3,
                                                                       33322, 16866, 33332, 4458,
                                                                       4467, 17220, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33942, 0, 3,
                                                                       33332, 16872, 33342, 4467,
                                                                       4476, 17238, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33972, 0, 3,
                                                                       33342, 16878, 33352, 4476,
                                                                       4485, 17256, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 34002, 0, 3,
                                                                       33352, 16884, 33362, 4485,
                                                                       4494, 17274, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34032, 0, 3,
                                                                       33372, 16896, 33402, 4512,
                                                                       4530, 17292, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34092, 0, 3,
                                                                       33402, 16914, 33432, 4530,
                                                                       4548, 17328, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34152, 0, 3,
                                                                       33432, 16932, 33462, 4548,
                                                                       4566, 17364, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34212, 0, 3,
                                                                       33462, 16950, 33492, 4566,
                                                                       4584, 17400, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34272, 0, 3,
                                                                       33492, 16968, 33522, 4584,
                                                                       4602, 17436, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34332, 0, 3,
                                                                       33522, 16986, 33552, 4602,
                                                                       4620, 17472, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34392, 0, 3,
                                                                       33552, 17004, 33582, 4620,
                                                                       4638, 17508, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34452, 0, 3,
                                                                       33582, 17022, 33612, 4638,
                                                                       4656, 17544, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34512, 0, 3,
                                                                       33612, 17040, 33642, 4656,
                                                                       4674, 17580, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34572, 0, 3,
                                                                       33642, 17058, 33672, 4674,
                                                                       4692, 17616, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34632, 0, 3,
                                                                       33702, 17094, 33732, 4728,
                                                                       4746, 17652, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34692, 0, 3,
                                                                       33732, 17112, 33762, 4746,
                                                                       4764, 17688, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34752, 0, 3,
                                                                       33762, 17130, 33792, 4764,
                                                                       4782, 17724, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34812, 0, 3,
                                                                       33792, 17148, 33822, 4782,
                                                                       4800, 17760, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34872, 0, 3,
                                                                       33822, 17166, 33852, 4800,
                                                                       4818, 17796, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34932, 0, 3,
                                                                       33852, 17184, 33882, 4818,
                                                                       4836, 17832, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34992, 0, 3,
                                                                       33882, 17202, 33912, 4836,
                                                                       4854, 17868, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 35052, 0, 3,
                                                                       33912, 17220, 33942, 4854,
                                                                       4872, 17904, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 35112, 0, 3,
                                                                       33942, 17238, 33972, 4872,
                                                                       4890, 17940, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 35172, 0, 3,
                                                                       33972, 17256, 34002, 4890,
                                                                       4908, 17976, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35232, 0, 3,
                                                                       34032, 17292, 34092, 4944,
                                                                       4974, 18012, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35332, 0, 3,
                                                                       34092, 17328, 34152, 4974,
                                                                       5004, 18072, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35432, 0, 3,
                                                                       34152, 17364, 34212, 5004,
                                                                       5034, 18132, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35532, 0, 3,
                                                                       34212, 17400, 34272, 5034,
                                                                       5064, 18192, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35632, 0, 3,
                                                                       34272, 17436, 34332, 5064,
                                                                       5094, 18252, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35732, 0, 3,
                                                                       34332, 17472, 34392, 5094,
                                                                       5124, 18312, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35832, 0, 3,
                                                                       34392, 17508, 34452, 5124,
                                                                       5154, 18372, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35932, 0, 3,
                                                                       34452, 17544, 34512, 5154,
                                                                       5184, 18432, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 36032, 0, 3,
                                                                       34512, 17580, 34572, 5184,
                                                                       5214, 18492, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 36132, 0, 3,
                                                                       34632, 17652, 34692, 5274,
                                                                       5304, 18552, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 36232, 0, 3,
                                                                       34692, 17688, 34752, 5304,
                                                                       5334, 18612, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 36332, 0, 3,
                                                                       34752, 17724, 34812, 5334,
                                                                       5364, 18672, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 36432, 0, 3,
                                                                       34812, 17760, 34872, 5364,
                                                                       5394, 18732, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 36532, 0, 3,
                                                                       34872, 17796, 34932, 5394,
                                                                       5424, 18792, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 36632, 0, 3,
                                                                       34932, 17832, 34992, 5424,
                                                                       5454, 18852, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 36732, 0, 3,
                                                                       34992, 17868, 35052, 5454,
                                                                       5484, 18912, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 36832, 0, 3,
                                                                       35052, 17904, 35112, 5484,
                                                                       5514, 18972, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 36932, 0, 3,
                                                                       35112, 17940, 35172, 5514,
                                                                       5544, 19032, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37032, 0, 3,
                                                                       35232, 18012, 35332, 5604,
                                                                       5649, 19092, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37182, 0, 3,
                                                                       35332, 18072, 35432, 5649,
                                                                       5694, 19182, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37332, 0, 3,
                                                                       35432, 18132, 35532, 5694,
                                                                       5739, 19272, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37482, 0, 3,
                                                                       35532, 18192, 35632, 5739,
                                                                       5784, 19362, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37632, 0, 3,
                                                                       35632, 18252, 35732, 5784,
                                                                       5829, 19452, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37782, 0, 3,
                                                                       35732, 18312, 35832, 5829,
                                                                       5874, 19542, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37932, 0, 3,
                                                                       35832, 18372, 35932, 5874,
                                                                       5919, 19632, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 38082, 0, 3,
                                                                       35932, 18432, 36032, 5919,
                                                                       5964, 19722, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 38232, 0, 3,
                                                                       36132, 18552, 36232, 6054,
                                                                       6099, 19812, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 38382, 0, 3,
                                                                       36232, 18612, 36332, 6099,
                                                                       6144, 19902, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 38532, 0, 3,
                                                                       36332, 18672, 36432, 6144,
                                                                       6189, 19992, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 38682, 0, 3,
                                                                       36432, 18732, 36532, 6189,
                                                                       6234, 20082, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 38832, 0, 3,
                                                                       36532, 18792, 36632, 6234,
                                                                       6279, 20172, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 38982, 0, 3,
                                                                       36632, 18852, 36732, 6279,
                                                                       6324, 20262, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 39132, 0, 3,
                                                                       36732, 18912, 36832, 6324,
                                                                       6369, 20352, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 39282, 0, 3,
                                                                       36832, 18972, 36932, 6369,
                                                                       6414, 20442, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39432, 0, 3,
                                                                       37032, 19092, 37182, 6504,
                                                                       6567, 20532, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39642, 0, 3,
                                                                       37182, 19182, 37332, 6567,
                                                                       6630, 20658, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39852, 0, 3,
                                                                       37332, 19272, 37482, 6630,
                                                                       6693, 20784, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 40062, 0, 3,
                                                                       37482, 19362, 37632, 6693,
                                                                       6756, 20910, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 40272, 0, 3,
                                                                       37632, 19452, 37782, 6756,
                                                                       6819, 21036, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 40482, 0, 3,
                                                                       37782, 19542, 37932, 6819,
                                                                       6882, 21162, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 40692, 0, 3,
                                                                       37932, 19632, 38082, 6882,
                                                                       6945, 21288, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 40902, 0, 3,
                                                                       38232, 19812, 38382, 7071,
                                                                       7134, 21414, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 41112, 0, 3,
                                                                       38382, 19902, 38532, 7134,
                                                                       7197, 21540, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 41322, 0, 3,
                                                                       38532, 19992, 38682, 7197,
                                                                       7260, 21666, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 41532, 0, 3,
                                                                       38682, 20082, 38832, 7260,
                                                                       7323, 21792, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 41742, 0, 3,
                                                                       38832, 20172, 38982, 7323,
                                                                       7386, 21918, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 41952, 0, 3,
                                                                       38982, 20262, 39132, 7386,
                                                                       7449, 22044, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 42162, 0, 3,
                                                                       39132, 20352, 39282, 7449,
                                                                       7512, 22170, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 42372, 0, 3,
                                                                       39432, 20532, 39642, 7638,
                                                                       7722, 22296, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 42652, 0, 3,
                                                                       39642, 20658, 39852, 7722,
                                                                       7806, 22464, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 42932, 0, 3,
                                                                       39852, 20784, 40062, 7806,
                                                                       7890, 22632, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 43212, 0, 3,
                                                                       40062, 20910, 40272, 7890,
                                                                       7974, 22800, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 43492, 0, 3,
                                                                       40272, 21036, 40482, 7974,
                                                                       8058, 22968, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 43772, 0, 3,
                                                                       40482, 21162, 40692, 8058,
                                                                       8142, 23136, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44052, 0, 3,
                                                                       40902, 21414, 41112, 8310,
                                                                       8394, 23304, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44332, 0, 3,
                                                                       41112, 21540, 41322, 8394,
                                                                       8478, 23472, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44612, 0, 3,
                                                                       41322, 21666, 41532, 8478,
                                                                       8562, 23640, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44892, 0, 3,
                                                                       41532, 21792, 41742, 8562,
                                                                       8646, 23808, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 45172, 0, 3,
                                                                       41742, 21918, 41952, 8646,
                                                                       8730, 23976, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 45452, 0, 3,
                                                                       41952, 22044, 42162, 8730,
                                                                       8814, 24144, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 45732, 0, 3,
                                                                       42372, 22296, 42652, 8982,
                                                                       9090, 24312, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46092, 0, 3,
                                                                       42652, 22464, 42932, 9090,
                                                                       9198, 24528, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46452, 0, 3,
                                                                       42932, 22632, 43212, 9198,
                                                                       9306, 24744, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46812, 0, 3,
                                                                       43212, 22800, 43492, 9306,
                                                                       9414, 24960, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47172, 0, 3,
                                                                       43492, 22968, 43772, 9414,
                                                                       9522, 25176, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47532, 0, 3,
                                                                       44052, 23304, 44332, 9738,
                                                                       9846, 25392, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47892, 0, 3,
                                                                       44332, 23472, 44612, 9846,
                                                                       9954, 25608, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 48252, 0, 3,
                                                                       44612, 23640, 44892, 9954,
                                                                       10062, 25824, ncols,
                                                                       gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 48612, 0, 3,
                                                                       44892, 23808, 45172,
                                                                       10062, 10170, 26040,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 48972, 0, 3,
                                                                       45172, 23976, 45452,
                                                                       10170, 10278, 26256,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 49332, 0, 3,
                                                                       45732, 24312, 46092,
                                                                       10494, 10629, 26472,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 49782, 0, 3,
                                                                       46092, 24528, 46452,
                                                                       10629, 10764, 26742,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 50232, 0, 3,
                                                                       46452, 24744, 46812,
                                                                       10764, 10899, 27012,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 50682, 0, 3,
                                                                       46812, 24960, 47172,
                                                                       10899, 11034, 27282,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 51132, 0, 3,
                                                                       47532, 25392, 47892,
                                                                       11304, 11439, 27552,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 51582, 0, 3,
                                                                       47892, 25608, 48252,
                                                                       11439, 11574, 27822,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 52032, 0, 3,
                                                                       48252, 25824, 48612,
                                                                       11574, 11709, 28092,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 52482, 0, 3,
                                                                       48612, 26040, 48972,
                                                                       11709, 11844, 28362,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 52932, 0, 3,
                                                                       49332, 26472, 49782,
                                                                       12114, 12279, 28632,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 53482, 0, 3,
                                                                       49782, 26742, 50232,
                                                                       12279, 12444, 28962,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 54032, 0, 3,
                                                                       50232, 27012, 50682,
                                                                       12444, 12609, 29292,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 54582, 0, 3,
                                                                       51132, 27552, 51582,
                                                                       12939, 13104, 29622,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 55132, 0, 3,
                                                                       51582, 27822, 52032,
                                                                       13104, 13269, 29952,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 55682, 0, 3,
                                                                       52032, 28092, 52482,
                                                                       13269, 13434, 30282,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 56232, 0, 3,
                                                                       52932, 28632, 53482,
                                                                       13764, 13962, 30612,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 56892, 0, 3,
                                                                       53482, 28962, 54032,
                                                                       13962, 14160, 31008,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 57552, 0, 3,
                                                                       54582, 29622, 55132,
                                                                       14556, 14754, 31404,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 58212, 0, 3,
                                                                       55132, 29952, 55682,
                                                                       14754, 14952, 31800,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 58872, 0, 3,
                                                                       56232, 30612, 56892,
                                                                       15348, 15582, 32196,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 59652, 0, 3,
                                                                       57552, 31404, 58212,
                                                                       16050, 16284, 32664,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 60432, 42372, 280, ncols);

                    simdfunc::contract_primitives(buffer, 60908, 44052, 280, ncols);

                    simdfunc::contract_primitives(buffer, 61384, 45732, 360, ncols);

                    simdfunc::contract_primitives(buffer, 61996, 47532, 360, ncols);

                    simdfunc::contract_primitives(buffer, 62608, 49332, 450, ncols);

                    simdfunc::contract_primitives(buffer, 63373, 51132, 450, ncols);

                    simdfunc::contract_primitives(buffer, 64138, 52932, 550, ncols);

                    simdfunc::contract_primitives(buffer, 65073, 54582, 550, ncols);

                    simdfunc::contract_primitives(buffer, 66008, 56232, 660, ncols);

                    simdfunc::contract_primitives(buffer, 67130, 57552, 660, ncols);

                    simdfunc::contract_primitives(buffer, 68252, 58872, 780, ncols);

                    simdfunc::contract_primitives(buffer, 69578, 59652, 780, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 60712, 60432, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 61188, 60908, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 61744, 61384, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 62356, 61996, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 63058, 62608, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 63823, 63373, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 64688, 64138, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 65623, 65073, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 66668, 66008, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 67790, 67130, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 69032, 68252, 78, 1, nmax);

        simdtrf::transform_f_inner(buffer, 70358, 69578, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 70904, 60712, 61744, 7, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 71492, 61188, 62356, 7, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 72080, 61744, 63058, 7, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 72836, 62356, 63823, 7, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 73592, 63058, 64688, 7, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 74537, 63823, 65623, 7, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 75482, 64688, 66668, 7, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 76637, 65623, 67790, 7, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 77792, 66668, 69032, 7, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 79178, 67790, 70358, 7, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 80564, 70904, 72080, 7, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 81740, 71492, 72836, 7, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 82916, 72080, 73592, 7, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 84428, 72836, 74537, 7, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 85940, 73592, 75482, 7, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 87830, 74537, 76637, 7, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 89720, 75482, 77792, 7, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 92030, 76637, 79178, 7, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 94340, 80564, 82916, 7, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 96300, 81740, 84428, 7, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 98260, 82916, 85940, 7, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 100780, 84428, 87830, 7,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 103300, 85940, 89720, 7,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 106450, 87830, 92030, 7,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 109600, 94340, 98260, 7,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 112540, 96300, 100780, 7,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 115480, 98260, 103300, 7,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 119260, 100780, 106450, 7,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 123040, 109600, 115480, 7,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 127156, 112540, 119260, 7,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 131272, 127156, 28, 7, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 131272, 77, nmax);

        simdtrf::transform_h_inner(buffer, 131272, 123040, 28, 7, nmax);

        simdtrf::transform_i_outer(values + 1001 * nvalues + n * npairs, nvalues, buffer, 131272,
                                   77, nmax);
    }

    for (size_t m = 0; m < 2002; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
