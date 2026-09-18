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


#include "SimdThreeCenterElectronRepulsionRsRecIHD.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
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
#include "SimdTransformD.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_ihd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ihd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 80200, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1430 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 80200, 28764, 6386, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 13,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 21, 3, 13,
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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4194, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4197, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4200, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4203, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4206, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4209, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4212, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4215, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4218, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4221, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4224, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4227, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4230, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4233, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4236, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4239, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4242, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4245, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4248, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4251, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4254, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4257, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4260, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4263, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4266, 3, 9, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4275, 3, 10, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4284, 3, 11, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4293, 3, 12, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4302, 3, 13, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4311, 3, 14, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4320, 3, 15, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4329, 3, 16, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4338, 3, 17, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4347, 3, 18, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4356, 3, 19, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4365, 3, 24, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4374, 3, 25, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4383, 3, 26, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4392, 3, 27, 90,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4401, 3, 28, 93,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4410, 3, 29, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4419, 3, 30, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4428, 3, 31, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4437, 3, 32, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4446, 3, 33, 108,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4455, 3, 34, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4464, 3, 42, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4482, 3, 45, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4500, 3, 48, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4518, 3, 51, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4536, 3, 54, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4554, 3, 57, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4572, 3, 60, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4590, 3, 63, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4608, 3, 66, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4626, 3, 69, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4644, 3, 81, 198,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4662, 3, 84, 204,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4680, 3, 87, 210,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4698, 3, 90, 216,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4716, 3, 93, 222,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4734, 3, 96, 228,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4752, 3, 99, 234,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4770, 3, 102, 240,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4788, 3, 105, 246,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4806, 3, 108, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4824, 3, 126, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4854, 3, 132, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4884, 3, 138, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4914, 3, 144, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4944, 3, 150, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4974, 3, 156, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5004, 3, 162, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5034, 3, 168, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5064, 3, 174, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5094, 3, 198, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5124, 3, 204, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5154, 3, 210, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5184, 3, 216, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5214, 3, 222, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5244, 3, 228, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5274, 3, 234, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5304, 3, 240, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5334, 3, 246, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5364, 3, 278, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5409, 3, 288, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5454, 3, 298, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5499, 3, 308, 553,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5544, 3, 318, 568,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5589, 3, 328, 583,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5634, 3, 338, 598,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5679, 3, 348, 613,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5724, 3, 388, 658,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5769, 3, 398, 673,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5814, 3, 408, 688,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5859, 3, 418, 703,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5904, 3, 428, 718,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5949, 3, 438, 733,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5994, 3, 448, 748,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6039, 3, 458, 763,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6084, 3, 508, 820,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6147, 3, 523, 841,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6210, 3, 538, 862,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6273, 3, 553, 883,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6336, 3, 568, 904,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6399, 3, 583, 925,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6462, 3, 598, 946,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6525, 3, 658,
                                                                       1009, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6588, 3, 673,
                                                                       1030, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6651, 3, 688,
                                                                       1051, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6714, 3, 703,
                                                                       1072, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6777, 3, 718,
                                                                       1093, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6840, 3, 733,
                                                                       1114, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6903, 3, 748,
                                                                       1135, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6966, 3, 820,
                                                                       1212, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7050, 3, 841,
                                                                       1240, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7134, 3, 862,
                                                                       1268, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7218, 3, 883,
                                                                       1296, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7302, 3, 904,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7386, 3, 925,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7470, 3, 1009,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7554, 3, 1030,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7638, 3, 1051,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7722, 3, 1072,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7806, 3, 1093,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7890, 3, 1114,
                                                                       1576, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7974, 3, 1212,
                                                                       1676, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8082, 3, 1240,
                                                                       1712, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8190, 3, 1268,
                                                                       1748, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8298, 3, 1296,
                                                                       1784, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8406, 3, 1324,
                                                                       1820, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8514, 3, 1436,
                                                                       1928, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8622, 3, 1464,
                                                                       1964, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8730, 3, 1492,
                                                                       2000, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8838, 3, 1520,
                                                                       2036, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8946, 3, 1548,
                                                                       2072, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9054, 3, 1676,
                                                                       2198, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9189, 3, 1712,
                                                                       2243, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9324, 3, 1748,
                                                                       2288, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9459, 3, 1784,
                                                                       2333, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9594, 3, 1928,
                                                                       2468, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9729, 3, 1964,
                                                                       2513, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9864, 3, 2000,
                                                                       2558, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9999, 3, 2036,
                                                                       2603, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10134, 3, 2198,
                                                                       2758, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10299, 3, 2243,
                                                                       2813, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10464, 3, 2288,
                                                                       2868, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10629, 3, 2468,
                                                                       3033, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10794, 3, 2513,
                                                                       3088, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10959, 3, 2558,
                                                                       3143, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11124, 3, 2758,
                                                                       3330, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11322, 3, 2813,
                                                                       3396, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11520, 3, 3033,
                                                                       3594, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11718, 3, 3088,
                                                                       3660, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 11916, 3, 3330,
                                                                       3882, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12150, 3, 3594,
                                                                       4116, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12384, 3, 7, 8,
                                                                       4194, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12390, 3, 8, 9,
                                                                       4197, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12396, 3, 9, 10,
                                                                       4200, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12402, 3, 10, 11,
                                                                       4203, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12408, 3, 11, 12,
                                                                       4206, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12414, 3, 12, 13,
                                                                       4209, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12420, 3, 13, 14,
                                                                       4212, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12426, 3, 14, 15,
                                                                       4215, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12432, 3, 15, 16,
                                                                       4218, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12438, 3, 16, 17,
                                                                       4221, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12444, 3, 17, 18,
                                                                       4224, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12450, 3, 18, 19,
                                                                       4227, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12456, 3, 22, 23,
                                                                       4230, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12462, 3, 23, 24,
                                                                       4233, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12468, 3, 24, 25,
                                                                       4236, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12474, 3, 25, 26,
                                                                       4239, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12480, 3, 26, 27,
                                                                       4242, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12486, 3, 27, 28,
                                                                       4245, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12492, 3, 28, 29,
                                                                       4248, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12498, 3, 29, 30,
                                                                       4251, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12504, 3, 30, 31,
                                                                       4254, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12510, 3, 31, 32,
                                                                       4257, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12516, 3, 32, 33,
                                                                       4260, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12522, 3, 33, 34,
                                                                       4263, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12528, 0, 3,
                                                                       12384, 4194, 12390, 4266,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12546, 0, 3,
                                                                       12390, 4197, 12396, 4275,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12564, 0, 3,
                                                                       12396, 4200, 12402, 4284,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12582, 0, 3,
                                                                       12402, 4203, 12408, 4293,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12600, 0, 3,
                                                                       12408, 4206, 12414, 4302,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12618, 0, 3,
                                                                       12414, 4209, 12420, 4311,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12636, 0, 3,
                                                                       12420, 4212, 12426, 4320,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12654, 0, 3,
                                                                       12426, 4215, 12432, 4329,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12672, 0, 3,
                                                                       12432, 4218, 12438, 4338,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12690, 0, 3,
                                                                       12438, 4221, 12444, 4347,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12708, 0, 3,
                                                                       12444, 4224, 12450, 4356,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12726, 0, 3,
                                                                       12456, 4230, 12462, 4365,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12744, 0, 3,
                                                                       12462, 4233, 12468, 4374,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12762, 0, 3,
                                                                       12468, 4236, 12474, 4383,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12780, 0, 3,
                                                                       12474, 4239, 12480, 4392,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12798, 0, 3,
                                                                       12480, 4242, 12486, 4401,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12816, 0, 3,
                                                                       12486, 4245, 12492, 4410,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12834, 0, 3,
                                                                       12492, 4248, 12498, 4419,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12852, 0, 3,
                                                                       12498, 4251, 12504, 4428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12870, 0, 3,
                                                                       12504, 4254, 12510, 4437,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12888, 0, 3,
                                                                       12510, 4257, 12516, 4446,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12906, 0, 3,
                                                                       12516, 4260, 12522, 4455,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12924, 0, 3,
                                                                       12528, 4266, 12546, 114,
                                                                       120, 4464, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12960, 0, 3,
                                                                       12546, 4275, 12564, 120,
                                                                       126, 4482, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12996, 0, 3,
                                                                       12564, 4284, 12582, 126,
                                                                       132, 4500, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13032, 0, 3,
                                                                       12582, 4293, 12600, 132,
                                                                       138, 4518, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13068, 0, 3,
                                                                       12600, 4302, 12618, 138,
                                                                       144, 4536, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13104, 0, 3,
                                                                       12618, 4311, 12636, 144,
                                                                       150, 4554, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13140, 0, 3,
                                                                       12636, 4320, 12654, 150,
                                                                       156, 4572, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13176, 0, 3,
                                                                       12654, 4329, 12672, 156,
                                                                       162, 4590, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13212, 0, 3,
                                                                       12672, 4338, 12690, 162,
                                                                       168, 4608, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13248, 0, 3,
                                                                       12690, 4347, 12708, 168,
                                                                       174, 4626, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13284, 0, 3,
                                                                       12726, 4365, 12744, 186,
                                                                       192, 4644, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13320, 0, 3,
                                                                       12744, 4374, 12762, 192,
                                                                       198, 4662, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13356, 0, 3,
                                                                       12762, 4383, 12780, 198,
                                                                       204, 4680, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13392, 0, 3,
                                                                       12780, 4392, 12798, 204,
                                                                       210, 4698, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13428, 0, 3,
                                                                       12798, 4401, 12816, 210,
                                                                       216, 4716, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13464, 0, 3,
                                                                       12816, 4410, 12834, 216,
                                                                       222, 4734, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13500, 0, 3,
                                                                       12834, 4419, 12852, 222,
                                                                       228, 4752, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13536, 0, 3,
                                                                       12852, 4428, 12870, 228,
                                                                       234, 4770, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13572, 0, 3,
                                                                       12870, 4437, 12888, 234,
                                                                       240, 4788, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13608, 0, 3,
                                                                       12888, 4446, 12906, 240,
                                                                       246, 4806, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13644, 0, 3,
                                                                       12924, 4464, 12960, 258,
                                                                       268, 4824, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13704, 0, 3,
                                                                       12960, 4482, 12996, 268,
                                                                       278, 4854, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13764, 0, 3,
                                                                       12996, 4500, 13032, 278,
                                                                       288, 4884, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13824, 0, 3,
                                                                       13032, 4518, 13068, 288,
                                                                       298, 4914, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13884, 0, 3,
                                                                       13068, 4536, 13104, 298,
                                                                       308, 4944, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13944, 0, 3,
                                                                       13104, 4554, 13140, 308,
                                                                       318, 4974, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14004, 0, 3,
                                                                       13140, 4572, 13176, 318,
                                                                       328, 5004, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14064, 0, 3,
                                                                       13176, 4590, 13212, 328,
                                                                       338, 5034, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14124, 0, 3,
                                                                       13212, 4608, 13248, 338,
                                                                       348, 5064, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14184, 0, 3,
                                                                       13284, 4644, 13320, 368,
                                                                       378, 5094, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14244, 0, 3,
                                                                       13320, 4662, 13356, 378,
                                                                       388, 5124, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14304, 0, 3,
                                                                       13356, 4680, 13392, 388,
                                                                       398, 5154, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14364, 0, 3,
                                                                       13392, 4698, 13428, 398,
                                                                       408, 5184, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14424, 0, 3,
                                                                       13428, 4716, 13464, 408,
                                                                       418, 5214, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14484, 0, 3,
                                                                       13464, 4734, 13500, 418,
                                                                       428, 5244, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14544, 0, 3,
                                                                       13500, 4752, 13536, 428,
                                                                       438, 5274, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14604, 0, 3,
                                                                       13536, 4770, 13572, 438,
                                                                       448, 5304, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14664, 0, 3,
                                                                       13572, 4788, 13608, 448,
                                                                       458, 5334, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14724, 0, 3,
                                                                       13644, 4824, 13704, 478,
                                                                       493, 5364, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14814, 0, 3,
                                                                       13704, 4854, 13764, 493,
                                                                       508, 5409, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14904, 0, 3,
                                                                       13764, 4884, 13824, 508,
                                                                       523, 5454, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14994, 0, 3,
                                                                       13824, 4914, 13884, 523,
                                                                       538, 5499, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15084, 0, 3,
                                                                       13884, 4944, 13944, 538,
                                                                       553, 5544, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15174, 0, 3,
                                                                       13944, 4974, 14004, 553,
                                                                       568, 5589, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15264, 0, 3,
                                                                       14004, 5004, 14064, 568,
                                                                       583, 5634, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15354, 0, 3,
                                                                       14064, 5034, 14124, 583,
                                                                       598, 5679, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15444, 0, 3,
                                                                       14184, 5094, 14244, 628,
                                                                       643, 5724, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15534, 0, 3,
                                                                       14244, 5124, 14304, 643,
                                                                       658, 5769, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15624, 0, 3,
                                                                       14304, 5154, 14364, 658,
                                                                       673, 5814, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15714, 0, 3,
                                                                       14364, 5184, 14424, 673,
                                                                       688, 5859, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15804, 0, 3,
                                                                       14424, 5214, 14484, 688,
                                                                       703, 5904, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15894, 0, 3,
                                                                       14484, 5244, 14544, 703,
                                                                       718, 5949, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15984, 0, 3,
                                                                       14544, 5274, 14604, 718,
                                                                       733, 5994, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16074, 0, 3,
                                                                       14604, 5304, 14664, 733,
                                                                       748, 6039, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16164, 0, 3,
                                                                       14724, 5364, 14814, 778,
                                                                       799, 6084, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16290, 0, 3,
                                                                       14814, 5409, 14904, 799,
                                                                       820, 6147, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16416, 0, 3,
                                                                       14904, 5454, 14994, 820,
                                                                       841, 6210, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16542, 0, 3,
                                                                       14994, 5499, 15084, 841,
                                                                       862, 6273, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16668, 0, 3,
                                                                       15084, 5544, 15174, 862,
                                                                       883, 6336, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16794, 0, 3,
                                                                       15174, 5589, 15264, 883,
                                                                       904, 6399, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16920, 0, 3,
                                                                       15264, 5634, 15354, 904,
                                                                       925, 6462, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17046, 0, 3,
                                                                       15444, 5724, 15534, 967,
                                                                       988, 6525, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17172, 0, 3,
                                                                       15534, 5769, 15624, 988,
                                                                       1009, 6588, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17298, 0, 3,
                                                                       15624, 5814, 15714, 1009,
                                                                       1030, 6651, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17424, 0, 3,
                                                                       15714, 5859, 15804, 1030,
                                                                       1051, 6714, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17550, 0, 3,
                                                                       15804, 5904, 15894, 1051,
                                                                       1072, 6777, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17676, 0, 3,
                                                                       15894, 5949, 15984, 1072,
                                                                       1093, 6840, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17802, 0, 3,
                                                                       15984, 5994, 16074, 1093,
                                                                       1114, 6903, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17928, 0, 3,
                                                                       16164, 6084, 16290, 1156,
                                                                       1184, 6966, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18096, 0, 3,
                                                                       16290, 6147, 16416, 1184,
                                                                       1212, 7050, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18264, 0, 3,
                                                                       16416, 6210, 16542, 1212,
                                                                       1240, 7134, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18432, 0, 3,
                                                                       16542, 6273, 16668, 1240,
                                                                       1268, 7218, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18600, 0, 3,
                                                                       16668, 6336, 16794, 1268,
                                                                       1296, 7302, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18768, 0, 3,
                                                                       16794, 6399, 16920, 1296,
                                                                       1324, 7386, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18936, 0, 3,
                                                                       17046, 6525, 17172, 1380,
                                                                       1408, 7470, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19104, 0, 3,
                                                                       17172, 6588, 17298, 1408,
                                                                       1436, 7554, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19272, 0, 3,
                                                                       17298, 6651, 17424, 1436,
                                                                       1464, 7638, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19440, 0, 3,
                                                                       17424, 6714, 17550, 1464,
                                                                       1492, 7722, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19608, 0, 3,
                                                                       17550, 6777, 17676, 1492,
                                                                       1520, 7806, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19776, 0, 3,
                                                                       17676, 6840, 17802, 1520,
                                                                       1548, 7890, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19944, 0, 3,
                                                                       17928, 6966, 18096, 1604,
                                                                       1640, 7974, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20160, 0, 3,
                                                                       18096, 7050, 18264, 1640,
                                                                       1676, 8082, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20376, 0, 3,
                                                                       18264, 7134, 18432, 1676,
                                                                       1712, 8190, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20592, 0, 3,
                                                                       18432, 7218, 18600, 1712,
                                                                       1748, 8298, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20808, 0, 3,
                                                                       18600, 7302, 18768, 1748,
                                                                       1784, 8406, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21024, 0, 3,
                                                                       18936, 7470, 19104, 1856,
                                                                       1892, 8514, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21240, 0, 3,
                                                                       19104, 7554, 19272, 1892,
                                                                       1928, 8622, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21456, 0, 3,
                                                                       19272, 7638, 19440, 1928,
                                                                       1964, 8730, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21672, 0, 3,
                                                                       19440, 7722, 19608, 1964,
                                                                       2000, 8838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21888, 0, 3,
                                                                       19608, 7806, 19776, 2000,
                                                                       2036, 8946, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22104, 0, 3,
                                                                       19944, 7974, 20160, 2108,
                                                                       2153, 9054, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22374, 0, 3,
                                                                       20160, 8082, 20376, 2153,
                                                                       2198, 9189, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22644, 0, 3,
                                                                       20376, 8190, 20592, 2198,
                                                                       2243, 9324, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22914, 0, 3,
                                                                       20592, 8298, 20808, 2243,
                                                                       2288, 9459, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23184, 0, 3,
                                                                       21024, 8514, 21240, 2378,
                                                                       2423, 9594, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23454, 0, 3,
                                                                       21240, 8622, 21456, 2423,
                                                                       2468, 9729, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23724, 0, 3,
                                                                       21456, 8730, 21672, 2468,
                                                                       2513, 9864, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23994, 0, 3,
                                                                       21672, 8838, 21888, 2513,
                                                                       2558, 9999, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24264, 0, 3,
                                                                       22104, 9054, 22374, 2648,
                                                                       2703, 10134, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24594, 0, 3,
                                                                       22374, 9189, 22644, 2703,
                                                                       2758, 10299, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24924, 0, 3,
                                                                       22644, 9324, 22914, 2758,
                                                                       2813, 10464, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25254, 0, 3,
                                                                       23184, 9594, 23454, 2923,
                                                                       2978, 10629, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25584, 0, 3,
                                                                       23454, 9729, 23724, 2978,
                                                                       3033, 10794, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25914, 0, 3,
                                                                       23724, 9864, 23994, 3033,
                                                                       3088, 10959, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26244, 0, 3,
                                                                       24264, 10134, 24594, 3198,
                                                                       3264, 11124, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26640, 0, 3,
                                                                       24594, 10299, 24924, 3264,
                                                                       3330, 11322, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27036, 0, 3,
                                                                       25254, 10629, 25584, 3462,
                                                                       3528, 11520, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27432, 0, 3,
                                                                       25584, 10794, 25914, 3528,
                                                                       3594, 11718, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 27828, 0, 3,
                                                                       26244, 11124, 26640, 3726,
                                                                       3804, 11916, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 28296, 0, 3,
                                                                       27036, 11520, 27432, 3960,
                                                                       4038, 12150, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 28764, 17928, 168, ncols);

                    simdfunc::contract_primitives(buffer, 29072, 18936, 168, ncols);

                    simdfunc::contract_primitives(buffer, 29380, 19944, 216, ncols);

                    simdfunc::contract_primitives(buffer, 29776, 21024, 216, ncols);

                    simdfunc::contract_primitives(buffer, 30172, 22104, 270, ncols);

                    simdfunc::contract_primitives(buffer, 30667, 23184, 270, ncols);

                    simdfunc::contract_primitives(buffer, 31162, 24264, 330, ncols);

                    simdfunc::contract_primitives(buffer, 31767, 25254, 330, ncols);

                    simdfunc::contract_primitives(buffer, 32372, 26244, 396, ncols);

                    simdfunc::contract_primitives(buffer, 33098, 27036, 396, ncols);

                    simdfunc::contract_primitives(buffer, 33824, 27828, 468, ncols);

                    simdfunc::contract_primitives(buffer, 34682, 28296, 468, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 28932, 28764, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 29240, 29072, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 29596, 29380, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 29992, 29776, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 30442, 30172, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 30937, 30667, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 31492, 31162, 55, 1, nmax);

        simdtrf::transform_d_inner(buffer, 32097, 31767, 55, 1, nmax);

        simdtrf::transform_d_inner(buffer, 32768, 32372, 66, 1, nmax);

        simdtrf::transform_d_inner(buffer, 33494, 33098, 66, 1, nmax);

        simdtrf::transform_d_inner(buffer, 34292, 33824, 78, 1, nmax);

        simdtrf::transform_d_inner(buffer, 35150, 34682, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 35540, 28932, 29596, 5, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 35960, 29240, 29992, 5, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 36380, 29596, 30442, 5, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 36920, 29992, 30937, 5, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 37460, 30442, 31492, 5, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 38135, 30937, 32097, 5, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 38810, 31492, 32768, 5, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 39635, 32097, 33494, 5, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 40460, 32768, 34292, 5, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 41450, 33494, 35150, 5, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 42440, 35540, 36380, 5, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 43280, 35960, 36920, 5, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 44120, 36380, 37460, 5, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 45200, 36920, 38135, 5, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 46280, 37460, 38810, 5, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 47630, 38135, 39635, 5, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 48980, 38810, 40460, 5, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 50630, 39635, 41450, 5, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 52280, 42440, 44120, 5, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 53680, 43280, 45200, 5, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 55080, 44120, 46280, 5, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 56880, 45200, 47630, 5, nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 58680, 46280, 48980, 5, nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 60930, 47630, 50630, 5, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 63180, 52280, 55080, 5, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 65280, 53680, 56880, 5, nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 67380, 55080, 58680, 5, nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 70080, 56880, 60930, 5, nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 72780, 63180, 67380, 5, nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 75720, 65280, 70080, 5, nmax);

        simdtrf::transform_h_inner(buffer, 78660, 75720, 28, 5, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 78660, 55, nmax);

        simdtrf::transform_h_inner(buffer, 78660, 72780, 28, 5, nmax);

        simdtrf::transform_i_outer(values + 715 * nvalues + n * npairs, nvalues, buffer, 78660,
                                   55, nmax);
    }

    for (size_t m = 0; m < 1430; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
