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


#include "SimdThreeCenterElectronRepulsionRsRecGFK.hpp"

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
#include "SimdTransferGD.hpp"
#include "SimdTransferGF.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferHD.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_gfk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_gfk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 146123, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1890 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 146123, 117608, 9660, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2108, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2111, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2114, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2117, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2120, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2123, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2126, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2129, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2132, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2135, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2138, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2141, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2144, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2147, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2150, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2153, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2156, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2159, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2162, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2165, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2168, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2171, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2174, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2177, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2180, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2183, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2186, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2189, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2192, 3, 7, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2201, 3, 8, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2210, 3, 9, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2219, 3, 10, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2228, 3, 11, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2237, 3, 12, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2246, 3, 13, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2255, 3, 14, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2264, 3, 15, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2273, 3, 16, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2282, 3, 17, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2291, 3, 18, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2300, 3, 19, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2309, 3, 22, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2318, 3, 23, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2327, 3, 24, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2336, 3, 25, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2345, 3, 26, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2354, 3, 27, 90,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2363, 3, 28, 93,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2372, 3, 29, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2381, 3, 30, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2390, 3, 31, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2399, 3, 32, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2408, 3, 33, 108,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2417, 3, 34, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2426, 3, 36, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2444, 3, 39, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2462, 3, 42, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2480, 3, 45, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2498, 3, 48, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2516, 3, 51, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2534, 3, 54, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2552, 3, 57, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2570, 3, 60, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2588, 3, 63, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2606, 3, 66, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2624, 3, 69, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2642, 3, 75, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2660, 3, 78, 192,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2678, 3, 81, 198,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2696, 3, 84, 204,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2714, 3, 87, 210,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2732, 3, 90, 216,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2750, 3, 93, 222,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2768, 3, 96, 228,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2786, 3, 99, 234,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2804, 3, 102, 240,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2822, 3, 105, 246,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2840, 3, 108, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2858, 3, 114, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2888, 3, 120, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2918, 3, 126, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2948, 3, 132, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2978, 3, 138, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3008, 3, 144, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3038, 3, 150, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3068, 3, 156, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3098, 3, 162, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3128, 3, 168, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3158, 3, 174, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3188, 3, 186, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3218, 3, 192, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3248, 3, 198, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3278, 3, 204, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3308, 3, 210, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3338, 3, 216, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3368, 3, 222, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3398, 3, 228, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3428, 3, 234, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3458, 3, 240, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3488, 3, 246, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3518, 3, 258, 478,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3563, 3, 268, 493,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3608, 3, 278, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3653, 3, 288, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3698, 3, 298, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3743, 3, 308, 553,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3788, 3, 318, 568,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3833, 3, 328, 583,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3878, 3, 338, 598,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3923, 3, 348, 613,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3968, 3, 368, 628,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4013, 3, 378, 643,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4058, 3, 388, 658,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4103, 3, 398, 673,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4148, 3, 408, 688,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4193, 3, 418, 703,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4238, 3, 428, 718,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4283, 3, 438, 733,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4328, 3, 448, 748,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4373, 3, 458, 763,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4418, 3, 478, 778,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4481, 3, 493, 799,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4544, 3, 508, 820,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4607, 3, 523, 841,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4670, 3, 538, 862,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4733, 3, 553, 883,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4796, 3, 568, 904,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4859, 3, 583, 925,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4922, 3, 598, 946,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4985, 3, 628, 967,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5048, 3, 643, 988,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5111, 3, 658,
                                                                       1009, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5174, 3, 673,
                                                                       1030, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5237, 3, 688,
                                                                       1051, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5300, 3, 703,
                                                                       1072, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5363, 3, 718,
                                                                       1093, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5426, 3, 733,
                                                                       1114, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5489, 3, 748,
                                                                       1135, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5552, 3, 778,
                                                                       1156, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5636, 3, 799,
                                                                       1184, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5720, 3, 820,
                                                                       1212, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5804, 3, 841,
                                                                       1240, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5888, 3, 862,
                                                                       1268, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5972, 3, 883,
                                                                       1296, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6056, 3, 904,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6140, 3, 925,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6224, 3, 967,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6308, 3, 988,
                                                                       1408, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6392, 3, 1009,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6476, 3, 1030,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6560, 3, 1051,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6644, 3, 1072,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6728, 3, 1093,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6812, 3, 1114,
                                                                       1576, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6896, 3, 1156,
                                                                       1604, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7004, 3, 1184,
                                                                       1640, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7112, 3, 1212,
                                                                       1676, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7220, 3, 1240,
                                                                       1712, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7328, 3, 1268,
                                                                       1748, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7436, 3, 1296,
                                                                       1784, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7544, 3, 1324,
                                                                       1820, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7652, 3, 1380,
                                                                       1856, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7760, 3, 1408,
                                                                       1892, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7868, 3, 1436,
                                                                       1928, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7976, 3, 1464,
                                                                       1964, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8084, 3, 1492,
                                                                       2000, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8192, 3, 1520,
                                                                       2036, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8300, 3, 1548,
                                                                       2072, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8408, 3, 7, 8,
                                                                       2114, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8414, 3, 8, 9,
                                                                       2117, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8420, 3, 9, 10,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8426, 3, 10, 11,
                                                                       2123, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8432, 3, 11, 12,
                                                                       2126, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8438, 3, 12, 13,
                                                                       2129, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8444, 3, 13, 14,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8450, 3, 14, 15,
                                                                       2135, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8456, 3, 15, 16,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8462, 3, 16, 17,
                                                                       2141, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8468, 3, 17, 18,
                                                                       2144, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8474, 3, 18, 19,
                                                                       2147, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8480, 3, 22, 23,
                                                                       2156, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8486, 3, 23, 24,
                                                                       2159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8492, 3, 24, 25,
                                                                       2162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8498, 3, 25, 26,
                                                                       2165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8504, 3, 26, 27,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8510, 3, 27, 28,
                                                                       2171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8516, 3, 28, 29,
                                                                       2174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8522, 3, 29, 30,
                                                                       2177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8528, 3, 30, 31,
                                                                       2180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8534, 3, 31, 32,
                                                                       2183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8540, 3, 32, 33,
                                                                       2186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8546, 3, 33, 34,
                                                                       2189, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8552, 0, 3, 8408,
                                                                       2114, 8414, 2210, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8570, 0, 3, 8414,
                                                                       2117, 8420, 2219, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8588, 0, 3, 8420,
                                                                       2120, 8426, 2228, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8606, 0, 3, 8426,
                                                                       2123, 8432, 2237, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8624, 0, 3, 8432,
                                                                       2126, 8438, 2246, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8642, 0, 3, 8438,
                                                                       2129, 8444, 2255, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8660, 0, 3, 8444,
                                                                       2132, 8450, 2264, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8678, 0, 3, 8450,
                                                                       2135, 8456, 2273, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8696, 0, 3, 8456,
                                                                       2138, 8462, 2282, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8714, 0, 3, 8462,
                                                                       2141, 8468, 2291, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8732, 0, 3, 8468,
                                                                       2144, 8474, 2300, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8750, 0, 3, 8480,
                                                                       2156, 8486, 2327, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8768, 0, 3, 8486,
                                                                       2159, 8492, 2336, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8786, 0, 3, 8492,
                                                                       2162, 8498, 2345, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8804, 0, 3, 8498,
                                                                       2165, 8504, 2354, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8822, 0, 3, 8504,
                                                                       2168, 8510, 2363, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8840, 0, 3, 8510,
                                                                       2171, 8516, 2372, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8858, 0, 3, 8516,
                                                                       2174, 8522, 2381, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8876, 0, 3, 8522,
                                                                       2177, 8528, 2390, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8894, 0, 3, 8528,
                                                                       2180, 8534, 2399, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8912, 0, 3, 8534,
                                                                       2183, 8540, 2408, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8930, 0, 3, 8540,
                                                                       2186, 8546, 2417, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8948, 0, 3, 8552,
                                                                       2210, 8570, 114, 120,
                                                                       2462, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8984, 0, 3, 8570,
                                                                       2219, 8588, 120, 126,
                                                                       2480, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9020, 0, 3, 8588,
                                                                       2228, 8606, 126, 132,
                                                                       2498, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9056, 0, 3, 8606,
                                                                       2237, 8624, 132, 138,
                                                                       2516, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9092, 0, 3, 8624,
                                                                       2246, 8642, 138, 144,
                                                                       2534, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9128, 0, 3, 8642,
                                                                       2255, 8660, 144, 150,
                                                                       2552, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9164, 0, 3, 8660,
                                                                       2264, 8678, 150, 156,
                                                                       2570, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9200, 0, 3, 8678,
                                                                       2273, 8696, 156, 162,
                                                                       2588, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9236, 0, 3, 8696,
                                                                       2282, 8714, 162, 168,
                                                                       2606, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9272, 0, 3, 8714,
                                                                       2291, 8732, 168, 174,
                                                                       2624, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9308, 0, 3, 8750,
                                                                       2327, 8768, 186, 192,
                                                                       2678, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9344, 0, 3, 8768,
                                                                       2336, 8786, 192, 198,
                                                                       2696, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9380, 0, 3, 8786,
                                                                       2345, 8804, 198, 204,
                                                                       2714, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9416, 0, 3, 8804,
                                                                       2354, 8822, 204, 210,
                                                                       2732, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9452, 0, 3, 8822,
                                                                       2363, 8840, 210, 216,
                                                                       2750, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9488, 0, 3, 8840,
                                                                       2372, 8858, 216, 222,
                                                                       2768, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9524, 0, 3, 8858,
                                                                       2381, 8876, 222, 228,
                                                                       2786, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9560, 0, 3, 8876,
                                                                       2390, 8894, 228, 234,
                                                                       2804, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9596, 0, 3, 8894,
                                                                       2399, 8912, 234, 240,
                                                                       2822, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9632, 0, 3, 8912,
                                                                       2408, 8930, 240, 246,
                                                                       2840, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9668, 0, 3, 8948,
                                                                       2462, 8984, 258, 268,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9728, 0, 3, 8984,
                                                                       2480, 9020, 268, 278,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9788, 0, 3, 9020,
                                                                       2498, 9056, 278, 288,
                                                                       2978, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9848, 0, 3, 9056,
                                                                       2516, 9092, 288, 298,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9908, 0, 3, 9092,
                                                                       2534, 9128, 298, 308,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9968, 0, 3, 9128,
                                                                       2552, 9164, 308, 318,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10028, 0, 3, 9164,
                                                                       2570, 9200, 318, 328,
                                                                       3098, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10088, 0, 3, 9200,
                                                                       2588, 9236, 328, 338,
                                                                       3128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10148, 0, 3, 9236,
                                                                       2606, 9272, 338, 348,
                                                                       3158, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10208, 0, 3, 9308,
                                                                       2678, 9344, 368, 378,
                                                                       3248, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10268, 0, 3, 9344,
                                                                       2696, 9380, 378, 388,
                                                                       3278, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10328, 0, 3, 9380,
                                                                       2714, 9416, 388, 398,
                                                                       3308, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10388, 0, 3, 9416,
                                                                       2732, 9452, 398, 408,
                                                                       3338, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10448, 0, 3, 9452,
                                                                       2750, 9488, 408, 418,
                                                                       3368, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10508, 0, 3, 9488,
                                                                       2768, 9524, 418, 428,
                                                                       3398, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10568, 0, 3, 9524,
                                                                       2786, 9560, 428, 438,
                                                                       3428, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10628, 0, 3, 9560,
                                                                       2804, 9596, 438, 448,
                                                                       3458, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10688, 0, 3, 9596,
                                                                       2822, 9632, 448, 458,
                                                                       3488, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10748, 0, 3, 9668,
                                                                       2918, 9728, 478, 493,
                                                                       3608, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10838, 0, 3, 9728,
                                                                       2948, 9788, 493, 508,
                                                                       3653, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10928, 0, 3, 9788,
                                                                       2978, 9848, 508, 523,
                                                                       3698, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11018, 0, 3, 9848,
                                                                       3008, 9908, 523, 538,
                                                                       3743, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11108, 0, 3, 9908,
                                                                       3038, 9968, 538, 553,
                                                                       3788, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11198, 0, 3, 9968,
                                                                       3068, 10028, 553, 568,
                                                                       3833, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11288, 0, 3,
                                                                       10028, 3098, 10088, 568,
                                                                       583, 3878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11378, 0, 3,
                                                                       10088, 3128, 10148, 583,
                                                                       598, 3923, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11468, 0, 3,
                                                                       10208, 3248, 10268, 628,
                                                                       643, 4058, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11558, 0, 3,
                                                                       10268, 3278, 10328, 643,
                                                                       658, 4103, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11648, 0, 3,
                                                                       10328, 3308, 10388, 658,
                                                                       673, 4148, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11738, 0, 3,
                                                                       10388, 3338, 10448, 673,
                                                                       688, 4193, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11828, 0, 3,
                                                                       10448, 3368, 10508, 688,
                                                                       703, 4238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11918, 0, 3,
                                                                       10508, 3398, 10568, 703,
                                                                       718, 4283, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12008, 0, 3,
                                                                       10568, 3428, 10628, 718,
                                                                       733, 4328, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12098, 0, 3,
                                                                       10628, 3458, 10688, 733,
                                                                       748, 4373, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12188, 0, 3,
                                                                       10748, 3608, 10838, 778,
                                                                       799, 4544, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12314, 0, 3,
                                                                       10838, 3653, 10928, 799,
                                                                       820, 4607, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12440, 0, 3,
                                                                       10928, 3698, 11018, 820,
                                                                       841, 4670, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12566, 0, 3,
                                                                       11018, 3743, 11108, 841,
                                                                       862, 4733, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12692, 0, 3,
                                                                       11108, 3788, 11198, 862,
                                                                       883, 4796, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12818, 0, 3,
                                                                       11198, 3833, 11288, 883,
                                                                       904, 4859, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12944, 0, 3,
                                                                       11288, 3878, 11378, 904,
                                                                       925, 4922, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13070, 0, 3,
                                                                       11468, 4058, 11558, 967,
                                                                       988, 5111, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13196, 0, 3,
                                                                       11558, 4103, 11648, 988,
                                                                       1009, 5174, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13322, 0, 3,
                                                                       11648, 4148, 11738, 1009,
                                                                       1030, 5237, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13448, 0, 3,
                                                                       11738, 4193, 11828, 1030,
                                                                       1051, 5300, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13574, 0, 3,
                                                                       11828, 4238, 11918, 1051,
                                                                       1072, 5363, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13700, 0, 3,
                                                                       11918, 4283, 12008, 1072,
                                                                       1093, 5426, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13826, 0, 3,
                                                                       12008, 4328, 12098, 1093,
                                                                       1114, 5489, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13952, 0, 3,
                                                                       12188, 4544, 12314, 1156,
                                                                       1184, 5720, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14120, 0, 3,
                                                                       12314, 4607, 12440, 1184,
                                                                       1212, 5804, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14288, 0, 3,
                                                                       12440, 4670, 12566, 1212,
                                                                       1240, 5888, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14456, 0, 3,
                                                                       12566, 4733, 12692, 1240,
                                                                       1268, 5972, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14624, 0, 3,
                                                                       12692, 4796, 12818, 1268,
                                                                       1296, 6056, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14792, 0, 3,
                                                                       12818, 4859, 12944, 1296,
                                                                       1324, 6140, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14960, 0, 3,
                                                                       13070, 5111, 13196, 1380,
                                                                       1408, 6392, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15128, 0, 3,
                                                                       13196, 5174, 13322, 1408,
                                                                       1436, 6476, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15296, 0, 3,
                                                                       13322, 5237, 13448, 1436,
                                                                       1464, 6560, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15464, 0, 3,
                                                                       13448, 5300, 13574, 1464,
                                                                       1492, 6644, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15632, 0, 3,
                                                                       13574, 5363, 13700, 1492,
                                                                       1520, 6728, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15800, 0, 3,
                                                                       13700, 5426, 13826, 1520,
                                                                       1548, 6812, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15968, 0, 3,
                                                                       13952, 5720, 14120, 1604,
                                                                       1640, 7112, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16184, 0, 3,
                                                                       14120, 5804, 14288, 1640,
                                                                       1676, 7220, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16400, 0, 3,
                                                                       14288, 5888, 14456, 1676,
                                                                       1712, 7328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16616, 0, 3,
                                                                       14456, 5972, 14624, 1712,
                                                                       1748, 7436, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16832, 0, 3,
                                                                       14624, 6056, 14792, 1748,
                                                                       1784, 7544, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17048, 0, 3,
                                                                       14960, 6392, 15128, 1856,
                                                                       1892, 7868, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17264, 0, 3,
                                                                       15128, 6476, 15296, 1892,
                                                                       1928, 7976, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17480, 0, 3,
                                                                       15296, 6560, 15464, 1928,
                                                                       1964, 8084, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17696, 0, 3,
                                                                       15464, 6644, 15632, 1964,
                                                                       2000, 8192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17912, 0, 3,
                                                                       15632, 6728, 15800, 2000,
                                                                       2036, 8300, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18128, 3, 2108,
                                                                       2111, 8408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18138, 3, 2111,
                                                                       2114, 8414, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18148, 3, 2114,
                                                                       2117, 8420, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18158, 3, 2117,
                                                                       2120, 8426, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18168, 3, 2120,
                                                                       2123, 8432, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18178, 3, 2123,
                                                                       2126, 8438, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18188, 3, 2126,
                                                                       2129, 8444, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18198, 3, 2129,
                                                                       2132, 8450, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18208, 3, 2132,
                                                                       2135, 8456, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18218, 3, 2135,
                                                                       2138, 8462, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18228, 3, 2138,
                                                                       2141, 8468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18238, 3, 2141,
                                                                       2144, 8474, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18248, 3, 2150,
                                                                       2153, 8480, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18258, 3, 2153,
                                                                       2156, 8486, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18268, 3, 2156,
                                                                       2159, 8492, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18278, 3, 2159,
                                                                       2162, 8498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18288, 3, 2162,
                                                                       2165, 8504, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18298, 3, 2165,
                                                                       2168, 8510, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18308, 3, 2168,
                                                                       2171, 8516, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18318, 3, 2171,
                                                                       2174, 8522, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18328, 3, 2174,
                                                                       2177, 8528, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18338, 3, 2177,
                                                                       2180, 8534, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18348, 3, 2180,
                                                                       2183, 8540, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18358, 3, 2183,
                                                                       2186, 8546, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18368, 0, 3,
                                                                       18128, 8408, 18138, 2192,
                                                                       2201, 8552, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18398, 0, 3,
                                                                       18138, 8414, 18148, 2201,
                                                                       2210, 8570, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18428, 0, 3,
                                                                       18148, 8420, 18158, 2210,
                                                                       2219, 8588, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18458, 0, 3,
                                                                       18158, 8426, 18168, 2219,
                                                                       2228, 8606, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18488, 0, 3,
                                                                       18168, 8432, 18178, 2228,
                                                                       2237, 8624, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18518, 0, 3,
                                                                       18178, 8438, 18188, 2237,
                                                                       2246, 8642, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18548, 0, 3,
                                                                       18188, 8444, 18198, 2246,
                                                                       2255, 8660, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18578, 0, 3,
                                                                       18198, 8450, 18208, 2255,
                                                                       2264, 8678, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18608, 0, 3,
                                                                       18208, 8456, 18218, 2264,
                                                                       2273, 8696, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18638, 0, 3,
                                                                       18218, 8462, 18228, 2273,
                                                                       2282, 8714, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18668, 0, 3,
                                                                       18228, 8468, 18238, 2282,
                                                                       2291, 8732, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18698, 0, 3,
                                                                       18248, 8480, 18258, 2309,
                                                                       2318, 8750, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18728, 0, 3,
                                                                       18258, 8486, 18268, 2318,
                                                                       2327, 8768, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18758, 0, 3,
                                                                       18268, 8492, 18278, 2327,
                                                                       2336, 8786, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18788, 0, 3,
                                                                       18278, 8498, 18288, 2336,
                                                                       2345, 8804, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18818, 0, 3,
                                                                       18288, 8504, 18298, 2345,
                                                                       2354, 8822, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18848, 0, 3,
                                                                       18298, 8510, 18308, 2354,
                                                                       2363, 8840, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18878, 0, 3,
                                                                       18308, 8516, 18318, 2363,
                                                                       2372, 8858, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18908, 0, 3,
                                                                       18318, 8522, 18328, 2372,
                                                                       2381, 8876, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18938, 0, 3,
                                                                       18328, 8528, 18338, 2381,
                                                                       2390, 8894, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18968, 0, 3,
                                                                       18338, 8534, 18348, 2390,
                                                                       2399, 8912, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18998, 0, 3,
                                                                       18348, 8540, 18358, 2399,
                                                                       2408, 8930, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19028, 0, 3,
                                                                       18368, 8552, 18398, 2426,
                                                                       2444, 8948, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19088, 0, 3,
                                                                       18398, 8570, 18428, 2444,
                                                                       2462, 8984, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19148, 0, 3,
                                                                       18428, 8588, 18458, 2462,
                                                                       2480, 9020, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19208, 0, 3,
                                                                       18458, 8606, 18488, 2480,
                                                                       2498, 9056, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19268, 0, 3,
                                                                       18488, 8624, 18518, 2498,
                                                                       2516, 9092, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19328, 0, 3,
                                                                       18518, 8642, 18548, 2516,
                                                                       2534, 9128, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19388, 0, 3,
                                                                       18548, 8660, 18578, 2534,
                                                                       2552, 9164, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19448, 0, 3,
                                                                       18578, 8678, 18608, 2552,
                                                                       2570, 9200, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19508, 0, 3,
                                                                       18608, 8696, 18638, 2570,
                                                                       2588, 9236, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19568, 0, 3,
                                                                       18638, 8714, 18668, 2588,
                                                                       2606, 9272, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19628, 0, 3,
                                                                       18698, 8750, 18728, 2642,
                                                                       2660, 9308, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19688, 0, 3,
                                                                       18728, 8768, 18758, 2660,
                                                                       2678, 9344, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19748, 0, 3,
                                                                       18758, 8786, 18788, 2678,
                                                                       2696, 9380, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19808, 0, 3,
                                                                       18788, 8804, 18818, 2696,
                                                                       2714, 9416, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19868, 0, 3,
                                                                       18818, 8822, 18848, 2714,
                                                                       2732, 9452, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19928, 0, 3,
                                                                       18848, 8840, 18878, 2732,
                                                                       2750, 9488, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19988, 0, 3,
                                                                       18878, 8858, 18908, 2750,
                                                                       2768, 9524, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20048, 0, 3,
                                                                       18908, 8876, 18938, 2768,
                                                                       2786, 9560, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20108, 0, 3,
                                                                       18938, 8894, 18968, 2786,
                                                                       2804, 9596, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20168, 0, 3,
                                                                       18968, 8912, 18998, 2804,
                                                                       2822, 9632, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20228, 0, 3,
                                                                       19028, 8948, 19088, 2858,
                                                                       2888, 9668, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20328, 0, 3,
                                                                       19088, 8984, 19148, 2888,
                                                                       2918, 9728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20428, 0, 3,
                                                                       19148, 9020, 19208, 2918,
                                                                       2948, 9788, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20528, 0, 3,
                                                                       19208, 9056, 19268, 2948,
                                                                       2978, 9848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20628, 0, 3,
                                                                       19268, 9092, 19328, 2978,
                                                                       3008, 9908, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20728, 0, 3,
                                                                       19328, 9128, 19388, 3008,
                                                                       3038, 9968, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20828, 0, 3,
                                                                       19388, 9164, 19448, 3038,
                                                                       3068, 10028, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20928, 0, 3,
                                                                       19448, 9200, 19508, 3068,
                                                                       3098, 10088, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21028, 0, 3,
                                                                       19508, 9236, 19568, 3098,
                                                                       3128, 10148, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21128, 0, 3,
                                                                       19628, 9308, 19688, 3188,
                                                                       3218, 10208, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21228, 0, 3,
                                                                       19688, 9344, 19748, 3218,
                                                                       3248, 10268, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21328, 0, 3,
                                                                       19748, 9380, 19808, 3248,
                                                                       3278, 10328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21428, 0, 3,
                                                                       19808, 9416, 19868, 3278,
                                                                       3308, 10388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21528, 0, 3,
                                                                       19868, 9452, 19928, 3308,
                                                                       3338, 10448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21628, 0, 3,
                                                                       19928, 9488, 19988, 3338,
                                                                       3368, 10508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21728, 0, 3,
                                                                       19988, 9524, 20048, 3368,
                                                                       3398, 10568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21828, 0, 3,
                                                                       20048, 9560, 20108, 3398,
                                                                       3428, 10628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21928, 0, 3,
                                                                       20108, 9596, 20168, 3428,
                                                                       3458, 10688, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22028, 0, 3,
                                                                       20228, 9668, 20328, 3518,
                                                                       3563, 10748, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22178, 0, 3,
                                                                       20328, 9728, 20428, 3563,
                                                                       3608, 10838, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22328, 0, 3,
                                                                       20428, 9788, 20528, 3608,
                                                                       3653, 10928, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22478, 0, 3,
                                                                       20528, 9848, 20628, 3653,
                                                                       3698, 11018, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22628, 0, 3,
                                                                       20628, 9908, 20728, 3698,
                                                                       3743, 11108, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22778, 0, 3,
                                                                       20728, 9968, 20828, 3743,
                                                                       3788, 11198, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22928, 0, 3,
                                                                       20828, 10028, 20928, 3788,
                                                                       3833, 11288, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23078, 0, 3,
                                                                       20928, 10088, 21028, 3833,
                                                                       3878, 11378, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23228, 0, 3,
                                                                       21128, 10208, 21228, 3968,
                                                                       4013, 11468, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23378, 0, 3,
                                                                       21228, 10268, 21328, 4013,
                                                                       4058, 11558, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23528, 0, 3,
                                                                       21328, 10328, 21428, 4058,
                                                                       4103, 11648, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23678, 0, 3,
                                                                       21428, 10388, 21528, 4103,
                                                                       4148, 11738, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23828, 0, 3,
                                                                       21528, 10448, 21628, 4148,
                                                                       4193, 11828, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23978, 0, 3,
                                                                       21628, 10508, 21728, 4193,
                                                                       4238, 11918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24128, 0, 3,
                                                                       21728, 10568, 21828, 4238,
                                                                       4283, 12008, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24278, 0, 3,
                                                                       21828, 10628, 21928, 4283,
                                                                       4328, 12098, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 24428, 0, 3,
                                                                       22028, 10748, 22178, 4418,
                                                                       4481, 12188, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 24638, 0, 3,
                                                                       22178, 10838, 22328, 4481,
                                                                       4544, 12314, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 24848, 0, 3,
                                                                       22328, 10928, 22478, 4544,
                                                                       4607, 12440, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25058, 0, 3,
                                                                       22478, 11018, 22628, 4607,
                                                                       4670, 12566, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25268, 0, 3,
                                                                       22628, 11108, 22778, 4670,
                                                                       4733, 12692, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25478, 0, 3,
                                                                       22778, 11198, 22928, 4733,
                                                                       4796, 12818, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25688, 0, 3,
                                                                       22928, 11288, 23078, 4796,
                                                                       4859, 12944, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25898, 0, 3,
                                                                       23228, 11468, 23378, 4985,
                                                                       5048, 13070, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26108, 0, 3,
                                                                       23378, 11558, 23528, 5048,
                                                                       5111, 13196, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26318, 0, 3,
                                                                       23528, 11648, 23678, 5111,
                                                                       5174, 13322, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26528, 0, 3,
                                                                       23678, 11738, 23828, 5174,
                                                                       5237, 13448, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26738, 0, 3,
                                                                       23828, 11828, 23978, 5237,
                                                                       5300, 13574, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26948, 0, 3,
                                                                       23978, 11918, 24128, 5300,
                                                                       5363, 13700, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 27158, 0, 3,
                                                                       24128, 12008, 24278, 5363,
                                                                       5426, 13826, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27368, 0, 3,
                                                                       24428, 12188, 24638, 5552,
                                                                       5636, 13952, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27648, 0, 3,
                                                                       24638, 12314, 24848, 5636,
                                                                       5720, 14120, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27928, 0, 3,
                                                                       24848, 12440, 25058, 5720,
                                                                       5804, 14288, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28208, 0, 3,
                                                                       25058, 12566, 25268, 5804,
                                                                       5888, 14456, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28488, 0, 3,
                                                                       25268, 12692, 25478, 5888,
                                                                       5972, 14624, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28768, 0, 3,
                                                                       25478, 12818, 25688, 5972,
                                                                       6056, 14792, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 29048, 0, 3,
                                                                       25898, 13070, 26108, 6224,
                                                                       6308, 14960, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 29328, 0, 3,
                                                                       26108, 13196, 26318, 6308,
                                                                       6392, 15128, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 29608, 0, 3,
                                                                       26318, 13322, 26528, 6392,
                                                                       6476, 15296, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 29888, 0, 3,
                                                                       26528, 13448, 26738, 6476,
                                                                       6560, 15464, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30168, 0, 3,
                                                                       26738, 13574, 26948, 6560,
                                                                       6644, 15632, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30448, 0, 3,
                                                                       26948, 13700, 27158, 6644,
                                                                       6728, 15800, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30728, 0, 3,
                                                                       27368, 13952, 27648, 6896,
                                                                       7004, 15968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 31088, 0, 3,
                                                                       27648, 14120, 27928, 7004,
                                                                       7112, 16184, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 31448, 0, 3,
                                                                       27928, 14288, 28208, 7112,
                                                                       7220, 16400, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 31808, 0, 3,
                                                                       28208, 14456, 28488, 7220,
                                                                       7328, 16616, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32168, 0, 3,
                                                                       28488, 14624, 28768, 7328,
                                                                       7436, 16832, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32528, 0, 3,
                                                                       29048, 14960, 29328, 7652,
                                                                       7760, 17048, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32888, 0, 3,
                                                                       29328, 15128, 29608, 7760,
                                                                       7868, 17264, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33248, 0, 3,
                                                                       29608, 15296, 29888, 7868,
                                                                       7976, 17480, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33608, 0, 3,
                                                                       29888, 15464, 30168, 7976,
                                                                       8084, 17696, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33968, 0, 3,
                                                                       30168, 15632, 30448, 8084,
                                                                       8192, 17912, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34328, 3, 8408,
                                                                       8414, 18148, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34343, 3, 8414,
                                                                       8420, 18158, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34358, 3, 8420,
                                                                       8426, 18168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34373, 3, 8426,
                                                                       8432, 18178, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34388, 3, 8432,
                                                                       8438, 18188, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34403, 3, 8438,
                                                                       8444, 18198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34418, 3, 8444,
                                                                       8450, 18208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34433, 3, 8450,
                                                                       8456, 18218, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34448, 3, 8456,
                                                                       8462, 18228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34463, 3, 8462,
                                                                       8468, 18238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34478, 3, 8480,
                                                                       8486, 18268, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34493, 3, 8486,
                                                                       8492, 18278, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34508, 3, 8492,
                                                                       8498, 18288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34523, 3, 8498,
                                                                       8504, 18298, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34538, 3, 8504,
                                                                       8510, 18308, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34553, 3, 8510,
                                                                       8516, 18318, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34568, 3, 8516,
                                                                       8522, 18328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34583, 3, 8522,
                                                                       8528, 18338, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34598, 3, 8528,
                                                                       8534, 18348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34613, 3, 8534,
                                                                       8540, 18358, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34628, 0, 3,
                                                                       34328, 18148, 34343, 8552,
                                                                       8570, 18428, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34673, 0, 3,
                                                                       34343, 18158, 34358, 8570,
                                                                       8588, 18458, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34718, 0, 3,
                                                                       34358, 18168, 34373, 8588,
                                                                       8606, 18488, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34763, 0, 3,
                                                                       34373, 18178, 34388, 8606,
                                                                       8624, 18518, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34808, 0, 3,
                                                                       34388, 18188, 34403, 8624,
                                                                       8642, 18548, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34853, 0, 3,
                                                                       34403, 18198, 34418, 8642,
                                                                       8660, 18578, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34898, 0, 3,
                                                                       34418, 18208, 34433, 8660,
                                                                       8678, 18608, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34943, 0, 3,
                                                                       34433, 18218, 34448, 8678,
                                                                       8696, 18638, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34988, 0, 3,
                                                                       34448, 18228, 34463, 8696,
                                                                       8714, 18668, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35033, 0, 3,
                                                                       34478, 18268, 34493, 8750,
                                                                       8768, 18758, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35078, 0, 3,
                                                                       34493, 18278, 34508, 8768,
                                                                       8786, 18788, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35123, 0, 3,
                                                                       34508, 18288, 34523, 8786,
                                                                       8804, 18818, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35168, 0, 3,
                                                                       34523, 18298, 34538, 8804,
                                                                       8822, 18848, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35213, 0, 3,
                                                                       34538, 18308, 34553, 8822,
                                                                       8840, 18878, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35258, 0, 3,
                                                                       34553, 18318, 34568, 8840,
                                                                       8858, 18908, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35303, 0, 3,
                                                                       34568, 18328, 34583, 8858,
                                                                       8876, 18938, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35348, 0, 3,
                                                                       34583, 18338, 34598, 8876,
                                                                       8894, 18968, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35393, 0, 3,
                                                                       34598, 18348, 34613, 8894,
                                                                       8912, 18998, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35438, 0, 3,
                                                                       34628, 18428, 34673, 8948,
                                                                       8984, 19148, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35528, 0, 3,
                                                                       34673, 18458, 34718, 8984,
                                                                       9020, 19208, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35618, 0, 3,
                                                                       34718, 18488, 34763, 9020,
                                                                       9056, 19268, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35708, 0, 3,
                                                                       34763, 18518, 34808, 9056,
                                                                       9092, 19328, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35798, 0, 3,
                                                                       34808, 18548, 34853, 9092,
                                                                       9128, 19388, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35888, 0, 3,
                                                                       34853, 18578, 34898, 9128,
                                                                       9164, 19448, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35978, 0, 3,
                                                                       34898, 18608, 34943, 9164,
                                                                       9200, 19508, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36068, 0, 3,
                                                                       34943, 18638, 34988, 9200,
                                                                       9236, 19568, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36158, 0, 3,
                                                                       35033, 18758, 35078, 9308,
                                                                       9344, 19748, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36248, 0, 3,
                                                                       35078, 18788, 35123, 9344,
                                                                       9380, 19808, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36338, 0, 3,
                                                                       35123, 18818, 35168, 9380,
                                                                       9416, 19868, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36428, 0, 3,
                                                                       35168, 18848, 35213, 9416,
                                                                       9452, 19928, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36518, 0, 3,
                                                                       35213, 18878, 35258, 9452,
                                                                       9488, 19988, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36608, 0, 3,
                                                                       35258, 18908, 35303, 9488,
                                                                       9524, 20048, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36698, 0, 3,
                                                                       35303, 18938, 35348, 9524,
                                                                       9560, 20108, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36788, 0, 3,
                                                                       35348, 18968, 35393, 9560,
                                                                       9596, 20168, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 36878, 0, 3,
                                                                       35438, 19148, 35528, 9668,
                                                                       9728, 20428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37028, 0, 3,
                                                                       35528, 19208, 35618, 9728,
                                                                       9788, 20528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37178, 0, 3,
                                                                       35618, 19268, 35708, 9788,
                                                                       9848, 20628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37328, 0, 3,
                                                                       35708, 19328, 35798, 9848,
                                                                       9908, 20728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37478, 0, 3,
                                                                       35798, 19388, 35888, 9908,
                                                                       9968, 20828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37628, 0, 3,
                                                                       35888, 19448, 35978, 9968,
                                                                       10028, 20928, ncols,
                                                                       gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37778, 0, 3,
                                                                       35978, 19508, 36068,
                                                                       10028, 10088, 21028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37928, 0, 3,
                                                                       36158, 19748, 36248,
                                                                       10208, 10268, 21328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38078, 0, 3,
                                                                       36248, 19808, 36338,
                                                                       10268, 10328, 21428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38228, 0, 3,
                                                                       36338, 19868, 36428,
                                                                       10328, 10388, 21528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38378, 0, 3,
                                                                       36428, 19928, 36518,
                                                                       10388, 10448, 21628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38528, 0, 3,
                                                                       36518, 19988, 36608,
                                                                       10448, 10508, 21728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38678, 0, 3,
                                                                       36608, 20048, 36698,
                                                                       10508, 10568, 21828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38828, 0, 3,
                                                                       36698, 20108, 36788,
                                                                       10568, 10628, 21928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 38978, 0, 3,
                                                                       36878, 20428, 37028,
                                                                       10748, 10838, 22328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 39203, 0, 3,
                                                                       37028, 20528, 37178,
                                                                       10838, 10928, 22478,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 39428, 0, 3,
                                                                       37178, 20628, 37328,
                                                                       10928, 11018, 22628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 39653, 0, 3,
                                                                       37328, 20728, 37478,
                                                                       11018, 11108, 22778,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 39878, 0, 3,
                                                                       37478, 20828, 37628,
                                                                       11108, 11198, 22928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40103, 0, 3,
                                                                       37628, 20928, 37778,
                                                                       11198, 11288, 23078,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40328, 0, 3,
                                                                       37928, 21328, 38078,
                                                                       11468, 11558, 23528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40553, 0, 3,
                                                                       38078, 21428, 38228,
                                                                       11558, 11648, 23678,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40778, 0, 3,
                                                                       38228, 21528, 38378,
                                                                       11648, 11738, 23828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41003, 0, 3,
                                                                       38378, 21628, 38528,
                                                                       11738, 11828, 23978,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41228, 0, 3,
                                                                       38528, 21728, 38678,
                                                                       11828, 11918, 24128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41453, 0, 3,
                                                                       38678, 21828, 38828,
                                                                       11918, 12008, 24278,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 41678, 0, 3,
                                                                       38978, 22328, 39203,
                                                                       12188, 12314, 24848,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 41993, 0, 3,
                                                                       39203, 22478, 39428,
                                                                       12314, 12440, 25058,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 42308, 0, 3,
                                                                       39428, 22628, 39653,
                                                                       12440, 12566, 25268,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 42623, 0, 3,
                                                                       39653, 22778, 39878,
                                                                       12566, 12692, 25478,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 42938, 0, 3,
                                                                       39878, 22928, 40103,
                                                                       12692, 12818, 25688,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 43253, 0, 3,
                                                                       40328, 23528, 40553,
                                                                       13070, 13196, 26318,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 43568, 0, 3,
                                                                       40553, 23678, 40778,
                                                                       13196, 13322, 26528,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 43883, 0, 3,
                                                                       40778, 23828, 41003,
                                                                       13322, 13448, 26738,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 44198, 0, 3,
                                                                       41003, 23978, 41228,
                                                                       13448, 13574, 26948,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 44513, 0, 3,
                                                                       41228, 24128, 41453,
                                                                       13574, 13700, 27158,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 44828, 0, 3,
                                                                       41678, 24848, 41993,
                                                                       13952, 14120, 27928,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 45248, 0, 3,
                                                                       41993, 25058, 42308,
                                                                       14120, 14288, 28208,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 45668, 0, 3,
                                                                       42308, 25268, 42623,
                                                                       14288, 14456, 28488,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 46088, 0, 3,
                                                                       42623, 25478, 42938,
                                                                       14456, 14624, 28768,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 46508, 0, 3,
                                                                       43253, 26318, 43568,
                                                                       14960, 15128, 29608,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 46928, 0, 3,
                                                                       43568, 26528, 43883,
                                                                       15128, 15296, 29888,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 47348, 0, 3,
                                                                       43883, 26738, 44198,
                                                                       15296, 15464, 30168,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 47768, 0, 3,
                                                                       44198, 26948, 44513,
                                                                       15464, 15632, 30448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 48188, 0, 3,
                                                                       44828, 27928, 45248,
                                                                       15968, 16184, 31448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 48728, 0, 3,
                                                                       45248, 28208, 45668,
                                                                       16184, 16400, 31808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 49268, 0, 3,
                                                                       45668, 28488, 46088,
                                                                       16400, 16616, 32168,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 49808, 0, 3,
                                                                       46508, 29608, 46928,
                                                                       17048, 17264, 33248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 50348, 0, 3,
                                                                       46928, 29888, 47348,
                                                                       17264, 17480, 33608,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 50888, 0, 3,
                                                                       47348, 30168, 47768,
                                                                       17480, 17696, 33968,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51428, 3, 18128,
                                                                       18138, 34328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51449, 3, 18138,
                                                                       18148, 34343, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51470, 3, 18148,
                                                                       18158, 34358, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51491, 3, 18158,
                                                                       18168, 34373, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51512, 3, 18168,
                                                                       18178, 34388, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51533, 3, 18178,
                                                                       18188, 34403, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51554, 3, 18188,
                                                                       18198, 34418, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51575, 3, 18198,
                                                                       18208, 34433, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51596, 3, 18208,
                                                                       18218, 34448, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51617, 3, 18218,
                                                                       18228, 34463, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51638, 3, 18248,
                                                                       18258, 34478, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51659, 3, 18258,
                                                                       18268, 34493, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51680, 3, 18268,
                                                                       18278, 34508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51701, 3, 18278,
                                                                       18288, 34523, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51722, 3, 18288,
                                                                       18298, 34538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51743, 3, 18298,
                                                                       18308, 34553, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51764, 3, 18308,
                                                                       18318, 34568, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51785, 3, 18318,
                                                                       18328, 34583, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51806, 3, 18328,
                                                                       18338, 34598, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51827, 3, 18338,
                                                                       18348, 34613, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51848, 0, 3,
                                                                       51428, 34328, 51449,
                                                                       18368, 18398, 34628,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51911, 0, 3,
                                                                       51449, 34343, 51470,
                                                                       18398, 18428, 34673,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51974, 0, 3,
                                                                       51470, 34358, 51491,
                                                                       18428, 18458, 34718,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52037, 0, 3,
                                                                       51491, 34373, 51512,
                                                                       18458, 18488, 34763,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52100, 0, 3,
                                                                       51512, 34388, 51533,
                                                                       18488, 18518, 34808,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52163, 0, 3,
                                                                       51533, 34403, 51554,
                                                                       18518, 18548, 34853,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52226, 0, 3,
                                                                       51554, 34418, 51575,
                                                                       18548, 18578, 34898,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52289, 0, 3,
                                                                       51575, 34433, 51596,
                                                                       18578, 18608, 34943,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52352, 0, 3,
                                                                       51596, 34448, 51617,
                                                                       18608, 18638, 34988,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52415, 0, 3,
                                                                       51638, 34478, 51659,
                                                                       18698, 18728, 35033,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52478, 0, 3,
                                                                       51659, 34493, 51680,
                                                                       18728, 18758, 35078,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52541, 0, 3,
                                                                       51680, 34508, 51701,
                                                                       18758, 18788, 35123,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52604, 0, 3,
                                                                       51701, 34523, 51722,
                                                                       18788, 18818, 35168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52667, 0, 3,
                                                                       51722, 34538, 51743,
                                                                       18818, 18848, 35213,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52730, 0, 3,
                                                                       51743, 34553, 51764,
                                                                       18848, 18878, 35258,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52793, 0, 3,
                                                                       51764, 34568, 51785,
                                                                       18878, 18908, 35303,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52856, 0, 3,
                                                                       51785, 34583, 51806,
                                                                       18908, 18938, 35348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52919, 0, 3,
                                                                       51806, 34598, 51827,
                                                                       18938, 18968, 35393,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52982, 0, 3,
                                                                       51848, 34628, 51911,
                                                                       19028, 19088, 35438,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53108, 0, 3,
                                                                       51911, 34673, 51974,
                                                                       19088, 19148, 35528,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53234, 0, 3,
                                                                       51974, 34718, 52037,
                                                                       19148, 19208, 35618,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53360, 0, 3,
                                                                       52037, 34763, 52100,
                                                                       19208, 19268, 35708,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53486, 0, 3,
                                                                       52100, 34808, 52163,
                                                                       19268, 19328, 35798,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53612, 0, 3,
                                                                       52163, 34853, 52226,
                                                                       19328, 19388, 35888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53738, 0, 3,
                                                                       52226, 34898, 52289,
                                                                       19388, 19448, 35978,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53864, 0, 3,
                                                                       52289, 34943, 52352,
                                                                       19448, 19508, 36068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53990, 0, 3,
                                                                       52415, 35033, 52478,
                                                                       19628, 19688, 36158,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 54116, 0, 3,
                                                                       52478, 35078, 52541,
                                                                       19688, 19748, 36248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 54242, 0, 3,
                                                                       52541, 35123, 52604,
                                                                       19748, 19808, 36338,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 54368, 0, 3,
                                                                       52604, 35168, 52667,
                                                                       19808, 19868, 36428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 54494, 0, 3,
                                                                       52667, 35213, 52730,
                                                                       19868, 19928, 36518,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 54620, 0, 3,
                                                                       52730, 35258, 52793,
                                                                       19928, 19988, 36608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 54746, 0, 3,
                                                                       52793, 35303, 52856,
                                                                       19988, 20048, 36698,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 54872, 0, 3,
                                                                       52856, 35348, 52919,
                                                                       20048, 20108, 36788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54998, 0, 3,
                                                                       52982, 35438, 53108,
                                                                       20228, 20328, 36878,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 55208, 0, 3,
                                                                       53108, 35528, 53234,
                                                                       20328, 20428, 37028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 55418, 0, 3,
                                                                       53234, 35618, 53360,
                                                                       20428, 20528, 37178,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 55628, 0, 3,
                                                                       53360, 35708, 53486,
                                                                       20528, 20628, 37328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 55838, 0, 3,
                                                                       53486, 35798, 53612,
                                                                       20628, 20728, 37478,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 56048, 0, 3,
                                                                       53612, 35888, 53738,
                                                                       20728, 20828, 37628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 56258, 0, 3,
                                                                       53738, 35978, 53864,
                                                                       20828, 20928, 37778,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 56468, 0, 3,
                                                                       53990, 36158, 54116,
                                                                       21128, 21228, 37928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 56678, 0, 3,
                                                                       54116, 36248, 54242,
                                                                       21228, 21328, 38078,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 56888, 0, 3,
                                                                       54242, 36338, 54368,
                                                                       21328, 21428, 38228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 57098, 0, 3,
                                                                       54368, 36428, 54494,
                                                                       21428, 21528, 38378,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 57308, 0, 3,
                                                                       54494, 36518, 54620,
                                                                       21528, 21628, 38528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 57518, 0, 3,
                                                                       54620, 36608, 54746,
                                                                       21628, 21728, 38678,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 57728, 0, 3,
                                                                       54746, 36698, 54872,
                                                                       21728, 21828, 38828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 57938, 0, 3,
                                                                       54998, 36878, 55208,
                                                                       22028, 22178, 38978,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 58253, 0, 3,
                                                                       55208, 37028, 55418,
                                                                       22178, 22328, 39203,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 58568, 0, 3,
                                                                       55418, 37178, 55628,
                                                                       22328, 22478, 39428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 58883, 0, 3,
                                                                       55628, 37328, 55838,
                                                                       22478, 22628, 39653,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 59198, 0, 3,
                                                                       55838, 37478, 56048,
                                                                       22628, 22778, 39878,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 59513, 0, 3,
                                                                       56048, 37628, 56258,
                                                                       22778, 22928, 40103,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 59828, 0, 3,
                                                                       56468, 37928, 56678,
                                                                       23228, 23378, 40328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 60143, 0, 3,
                                                                       56678, 38078, 56888,
                                                                       23378, 23528, 40553,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 60458, 0, 3,
                                                                       56888, 38228, 57098,
                                                                       23528, 23678, 40778,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 60773, 0, 3,
                                                                       57098, 38378, 57308,
                                                                       23678, 23828, 41003,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 61088, 0, 3,
                                                                       57308, 38528, 57518,
                                                                       23828, 23978, 41228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 61403, 0, 3,
                                                                       57518, 38678, 57728,
                                                                       23978, 24128, 41453,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 61718, 0, 3,
                                                                       57938, 38978, 58253,
                                                                       24428, 24638, 41678,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 62159, 0, 3,
                                                                       58253, 39203, 58568,
                                                                       24638, 24848, 41993,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 62600, 0, 3,
                                                                       58568, 39428, 58883,
                                                                       24848, 25058, 42308,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 63041, 0, 3,
                                                                       58883, 39653, 59198,
                                                                       25058, 25268, 42623,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 63482, 0, 3,
                                                                       59198, 39878, 59513,
                                                                       25268, 25478, 42938,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 63923, 0, 3,
                                                                       59828, 40328, 60143,
                                                                       25898, 26108, 43253,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 64364, 0, 3,
                                                                       60143, 40553, 60458,
                                                                       26108, 26318, 43568,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 64805, 0, 3,
                                                                       60458, 40778, 60773,
                                                                       26318, 26528, 43883,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 65246, 0, 3,
                                                                       60773, 41003, 61088,
                                                                       26528, 26738, 44198,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 65687, 0, 3,
                                                                       61088, 41228, 61403,
                                                                       26738, 26948, 44513,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 66128, 0, 3,
                                                                       61718, 41678, 62159,
                                                                       27368, 27648, 44828,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 66716, 0, 3,
                                                                       62159, 41993, 62600,
                                                                       27648, 27928, 45248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 67304, 0, 3,
                                                                       62600, 42308, 63041,
                                                                       27928, 28208, 45668,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 67892, 0, 3,
                                                                       63041, 42623, 63482,
                                                                       28208, 28488, 46088,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 68480, 0, 3,
                                                                       63923, 43253, 64364,
                                                                       29048, 29328, 46508,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 69068, 0, 3,
                                                                       64364, 43568, 64805,
                                                                       29328, 29608, 46928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 69656, 0, 3,
                                                                       64805, 43883, 65246,
                                                                       29608, 29888, 47348,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 70244, 0, 3,
                                                                       65246, 44198, 65687,
                                                                       29888, 30168, 47768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 70832, 0, 3,
                                                                       66128, 44828, 66716,
                                                                       30728, 31088, 48188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 71588, 0, 3,
                                                                       66716, 45248, 67304,
                                                                       31088, 31448, 48728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 72344, 0, 3,
                                                                       67304, 45668, 67892,
                                                                       31448, 31808, 49268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 73100, 0, 3,
                                                                       68480, 46508, 69068,
                                                                       32528, 32888, 49808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 73856, 0, 3,
                                                                       69068, 46928, 69656,
                                                                       32888, 33248, 50348,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 74612, 0, 3,
                                                                       69656, 47348, 70244,
                                                                       33248, 33608, 50888,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75368, 3, 34328,
                                                                       34343, 51470, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75396, 3, 34343,
                                                                       34358, 51491, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75424, 3, 34358,
                                                                       34373, 51512, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75452, 3, 34373,
                                                                       34388, 51533, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75480, 3, 34388,
                                                                       34403, 51554, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75508, 3, 34403,
                                                                       34418, 51575, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75536, 3, 34418,
                                                                       34433, 51596, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75564, 3, 34433,
                                                                       34448, 51617, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75592, 3, 34478,
                                                                       34493, 51680, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75620, 3, 34493,
                                                                       34508, 51701, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75648, 3, 34508,
                                                                       34523, 51722, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75676, 3, 34523,
                                                                       34538, 51743, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75704, 3, 34538,
                                                                       34553, 51764, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75732, 3, 34553,
                                                                       34568, 51785, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75760, 3, 34568,
                                                                       34583, 51806, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75788, 3, 34583,
                                                                       34598, 51827, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75816, 0, 3,
                                                                       75368, 51470, 75396,
                                                                       34628, 34673, 51974,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75900, 0, 3,
                                                                       75396, 51491, 75424,
                                                                       34673, 34718, 52037,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75984, 0, 3,
                                                                       75424, 51512, 75452,
                                                                       34718, 34763, 52100,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76068, 0, 3,
                                                                       75452, 51533, 75480,
                                                                       34763, 34808, 52163,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76152, 0, 3,
                                                                       75480, 51554, 75508,
                                                                       34808, 34853, 52226,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76236, 0, 3,
                                                                       75508, 51575, 75536,
                                                                       34853, 34898, 52289,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76320, 0, 3,
                                                                       75536, 51596, 75564,
                                                                       34898, 34943, 52352,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76404, 0, 3,
                                                                       75592, 51680, 75620,
                                                                       35033, 35078, 52541,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76488, 0, 3,
                                                                       75620, 51701, 75648,
                                                                       35078, 35123, 52604,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76572, 0, 3,
                                                                       75648, 51722, 75676,
                                                                       35123, 35168, 52667,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76656, 0, 3,
                                                                       75676, 51743, 75704,
                                                                       35168, 35213, 52730,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76740, 0, 3,
                                                                       75704, 51764, 75732,
                                                                       35213, 35258, 52793,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76824, 0, 3,
                                                                       75732, 51785, 75760,
                                                                       35258, 35303, 52856,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76908, 0, 3,
                                                                       75760, 51806, 75788,
                                                                       35303, 35348, 52919,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 76992, 0, 3,
                                                                       75816, 51974, 75900,
                                                                       35438, 35528, 53234,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 77160, 0, 3,
                                                                       75900, 52037, 75984,
                                                                       35528, 35618, 53360,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 77328, 0, 3,
                                                                       75984, 52100, 76068,
                                                                       35618, 35708, 53486,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 77496, 0, 3,
                                                                       76068, 52163, 76152,
                                                                       35708, 35798, 53612,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 77664, 0, 3,
                                                                       76152, 52226, 76236,
                                                                       35798, 35888, 53738,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 77832, 0, 3,
                                                                       76236, 52289, 76320,
                                                                       35888, 35978, 53864,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 78000, 0, 3,
                                                                       76404, 52541, 76488,
                                                                       36158, 36248, 54242,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 78168, 0, 3,
                                                                       76488, 52604, 76572,
                                                                       36248, 36338, 54368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 78336, 0, 3,
                                                                       76572, 52667, 76656,
                                                                       36338, 36428, 54494,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 78504, 0, 3,
                                                                       76656, 52730, 76740,
                                                                       36428, 36518, 54620,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 78672, 0, 3,
                                                                       76740, 52793, 76824,
                                                                       36518, 36608, 54746,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 78840, 0, 3,
                                                                       76824, 52856, 76908,
                                                                       36608, 36698, 54872,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 79008, 0, 3,
                                                                       76992, 53234, 77160,
                                                                       36878, 37028, 55418,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 79288, 0, 3,
                                                                       77160, 53360, 77328,
                                                                       37028, 37178, 55628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 79568, 0, 3,
                                                                       77328, 53486, 77496,
                                                                       37178, 37328, 55838,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 79848, 0, 3,
                                                                       77496, 53612, 77664,
                                                                       37328, 37478, 56048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 80128, 0, 3,
                                                                       77664, 53738, 77832,
                                                                       37478, 37628, 56258,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 80408, 0, 3,
                                                                       78000, 54242, 78168,
                                                                       37928, 38078, 56888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 80688, 0, 3,
                                                                       78168, 54368, 78336,
                                                                       38078, 38228, 57098,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 80968, 0, 3,
                                                                       78336, 54494, 78504,
                                                                       38228, 38378, 57308,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 81248, 0, 3,
                                                                       78504, 54620, 78672,
                                                                       38378, 38528, 57518,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 81528, 0, 3,
                                                                       78672, 54746, 78840,
                                                                       38528, 38678, 57728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 81808, 0, 3,
                                                                       79008, 55418, 79288,
                                                                       38978, 39203, 58568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 82228, 0, 3,
                                                                       79288, 55628, 79568,
                                                                       39203, 39428, 58883,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 82648, 0, 3,
                                                                       79568, 55838, 79848,
                                                                       39428, 39653, 59198,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 83068, 0, 3,
                                                                       79848, 56048, 80128,
                                                                       39653, 39878, 59513,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 83488, 0, 3,
                                                                       80408, 56888, 80688,
                                                                       40328, 40553, 60458,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 83908, 0, 3,
                                                                       80688, 57098, 80968,
                                                                       40553, 40778, 60773,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 84328, 0, 3,
                                                                       80968, 57308, 81248,
                                                                       40778, 41003, 61088,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 84748, 0, 3,
                                                                       81248, 57518, 81528,
                                                                       41003, 41228, 61403,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 85168, 0, 3,
                                                                       81808, 58568, 82228,
                                                                       41678, 41993, 62600,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 85756, 0, 3,
                                                                       82228, 58883, 82648,
                                                                       41993, 42308, 63041,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 86344, 0, 3,
                                                                       82648, 59198, 83068,
                                                                       42308, 42623, 63482,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 86932, 0, 3,
                                                                       83488, 60458, 83908,
                                                                       43253, 43568, 64805,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 87520, 0, 3,
                                                                       83908, 60773, 84328,
                                                                       43568, 43883, 65246,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 88108, 0, 3,
                                                                       84328, 61088, 84748,
                                                                       43883, 44198, 65687,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 88696, 0, 3,
                                                                       85168, 62600, 85756,
                                                                       44828, 45248, 67304,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 89480, 0, 3,
                                                                       85756, 63041, 86344,
                                                                       45248, 45668, 67892,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 90264, 0, 3,
                                                                       86932, 64805, 87520,
                                                                       46508, 46928, 69656,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 91048, 0, 3,
                                                                       87520, 65246, 88108,
                                                                       46928, 47348, 70244,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 91832, 0, 3,
                                                                       88696, 67304, 89480,
                                                                       48188, 48728, 72344,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 92840, 0, 3,
                                                                       90264, 69656, 91048,
                                                                       49808, 50348, 74612,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 93848, 3, 51428,
                                                                       51449, 75368, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 93884, 3, 51449,
                                                                       51470, 75396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 93920, 3, 51470,
                                                                       51491, 75424, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 93956, 3, 51491,
                                                                       51512, 75452, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 93992, 3, 51512,
                                                                       51533, 75480, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94028, 3, 51533,
                                                                       51554, 75508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94064, 3, 51554,
                                                                       51575, 75536, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94100, 3, 51575,
                                                                       51596, 75564, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94136, 3, 51638,
                                                                       51659, 75592, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94172, 3, 51659,
                                                                       51680, 75620, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94208, 3, 51680,
                                                                       51701, 75648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94244, 3, 51701,
                                                                       51722, 75676, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94280, 3, 51722,
                                                                       51743, 75704, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94316, 3, 51743,
                                                                       51764, 75732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94352, 3, 51764,
                                                                       51785, 75760, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94388, 3, 51785,
                                                                       51806, 75788, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 94424, 0, 3,
                                                                       93848, 75368, 93884,
                                                                       51848, 51911, 75816,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 94532, 0, 3,
                                                                       93884, 75396, 93920,
                                                                       51911, 51974, 75900,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 94640, 0, 3,
                                                                       93920, 75424, 93956,
                                                                       51974, 52037, 75984,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 94748, 0, 3,
                                                                       93956, 75452, 93992,
                                                                       52037, 52100, 76068,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 94856, 0, 3,
                                                                       93992, 75480, 94028,
                                                                       52100, 52163, 76152,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 94964, 0, 3,
                                                                       94028, 75508, 94064,
                                                                       52163, 52226, 76236,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95072, 0, 3,
                                                                       94064, 75536, 94100,
                                                                       52226, 52289, 76320,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95180, 0, 3,
                                                                       94136, 75592, 94172,
                                                                       52415, 52478, 76404,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95288, 0, 3,
                                                                       94172, 75620, 94208,
                                                                       52478, 52541, 76488,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95396, 0, 3,
                                                                       94208, 75648, 94244,
                                                                       52541, 52604, 76572,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95504, 0, 3,
                                                                       94244, 75676, 94280,
                                                                       52604, 52667, 76656,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95612, 0, 3,
                                                                       94280, 75704, 94316,
                                                                       52667, 52730, 76740,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95720, 0, 3,
                                                                       94316, 75732, 94352,
                                                                       52730, 52793, 76824,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95828, 0, 3,
                                                                       94352, 75760, 94388,
                                                                       52793, 52856, 76908,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 95936, 0, 3,
                                                                       94424, 75816, 94532,
                                                                       52982, 53108, 76992,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 96152, 0, 3,
                                                                       94532, 75900, 94640,
                                                                       53108, 53234, 77160,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 96368, 0, 3,
                                                                       94640, 75984, 94748,
                                                                       53234, 53360, 77328,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 96584, 0, 3,
                                                                       94748, 76068, 94856,
                                                                       53360, 53486, 77496,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 96800, 0, 3,
                                                                       94856, 76152, 94964,
                                                                       53486, 53612, 77664,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97016, 0, 3,
                                                                       94964, 76236, 95072,
                                                                       53612, 53738, 77832,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97232, 0, 3,
                                                                       95180, 76404, 95288,
                                                                       53990, 54116, 78000,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97448, 0, 3,
                                                                       95288, 76488, 95396,
                                                                       54116, 54242, 78168,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97664, 0, 3,
                                                                       95396, 76572, 95504,
                                                                       54242, 54368, 78336,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97880, 0, 3,
                                                                       95504, 76656, 95612,
                                                                       54368, 54494, 78504,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 98096, 0, 3,
                                                                       95612, 76740, 95720,
                                                                       54494, 54620, 78672,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 98312, 0, 3,
                                                                       95720, 76824, 95828,
                                                                       54620, 54746, 78840,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 98528, 0, 3,
                                                                       95936, 76992, 96152,
                                                                       54998, 55208, 79008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 98888, 0, 3,
                                                                       96152, 77160, 96368,
                                                                       55208, 55418, 79288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 99248, 0, 3,
                                                                       96368, 77328, 96584,
                                                                       55418, 55628, 79568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 99608, 0, 3,
                                                                       96584, 77496, 96800,
                                                                       55628, 55838, 79848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 99968, 0, 3,
                                                                       96800, 77664, 97016,
                                                                       55838, 56048, 80128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 100328, 0, 3,
                                                                       97232, 78000, 97448,
                                                                       56468, 56678, 80408,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 100688, 0, 3,
                                                                       97448, 78168, 97664,
                                                                       56678, 56888, 80688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 101048, 0, 3,
                                                                       97664, 78336, 97880,
                                                                       56888, 57098, 80968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 101408, 0, 3,
                                                                       97880, 78504, 98096,
                                                                       57098, 57308, 81248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 101768, 0, 3,
                                                                       98096, 78672, 98312,
                                                                       57308, 57518, 81528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 102128, 0, 3,
                                                                       98528, 79008, 98888,
                                                                       57938, 58253, 81808,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 102668, 0, 3,
                                                                       98888, 79288, 99248,
                                                                       58253, 58568, 82228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 103208, 0, 3,
                                                                       99248, 79568, 99608,
                                                                       58568, 58883, 82648,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 103748, 0, 3,
                                                                       99608, 79848, 99968,
                                                                       58883, 59198, 83068,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 104288, 0, 3,
                                                                       100328, 80408, 100688,
                                                                       59828, 60143, 83488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 104828, 0, 3,
                                                                       100688, 80688, 101048,
                                                                       60143, 60458, 83908,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 105368, 0, 3,
                                                                       101048, 80968, 101408,
                                                                       60458, 60773, 84328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 105908, 0, 3,
                                                                       101408, 81248, 101768,
                                                                       60773, 61088, 84748,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 106448, 0, 3,
                                                                       102128, 81808, 102668,
                                                                       61718, 62159, 85168,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 107204, 0, 3,
                                                                       102668, 82228, 103208,
                                                                       62159, 62600, 85756,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 107960, 0, 3,
                                                                       103208, 82648, 103748,
                                                                       62600, 63041, 86344,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 108716, 0, 3,
                                                                       104288, 83488, 104828,
                                                                       63923, 64364, 86932,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 109472, 0, 3,
                                                                       104828, 83908, 105368,
                                                                       64364, 64805, 87520,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 110228, 0, 3,
                                                                       105368, 84328, 105908,
                                                                       64805, 65246, 88108,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 110984, 0, 3,
                                                                       106448, 85168, 107204,
                                                                       66128, 66716, 88696,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 111992, 0, 3,
                                                                       107204, 85756, 107960,
                                                                       66716, 67304, 89480,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 113000, 0, 3,
                                                                       108716, 86932, 109472,
                                                                       68480, 69068, 90264,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 114008, 0, 3,
                                                                       109472, 87520, 110228,
                                                                       69068, 69656, 91048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 115016, 0, 3,
                                                                       110984, 88696, 111992,
                                                                       70832, 71588, 91832,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 116312, 0, 3,
                                                                       113000, 90264, 114008,
                                                                       73100, 73856, 92840,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 117608, 102128, 540, ncols);

                    simdfunc::contract_primitives(buffer, 118373, 104288, 540, ncols);

                    simdfunc::contract_primitives(buffer, 119138, 106448, 756, ncols);

                    simdfunc::contract_primitives(buffer, 120209, 108716, 756, ncols);

                    simdfunc::contract_primitives(buffer, 121280, 110984, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 122708, 113000, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 124136, 115016, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 125972, 116312, 1296, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 118148, 117608, 15, 1, nmax);

        simdtrf::transform_k_inner(buffer, 118913, 118373, 15, 1, nmax);

        simdtrf::transform_k_inner(buffer, 119894, 119138, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 120965, 120209, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 122288, 121280, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 123716, 122708, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 125432, 124136, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 127268, 125972, 36, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 127808, 118148, 119894, 15,
                                             nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 128483, 118913, 120965, 15,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 129158, 119894, 122288, 15,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 130103, 120965, 123716, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 131048, 122288, 125432, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 132308, 123716, 127268, 15,
                                             nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 133568, 127808, 129158, 15,
                                             nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 134918, 128483, 130103, 15,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 136268, 129158, 131048, 15,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 138158, 130103, 132308, 15,
                                             nmax);

        simdtrf::compute_hrr_gf_out_of_first(buffer, coordinates, 140048, 133568, 136268, 15,
                                             nmax);

        simdtrf::compute_hrr_gf_out_of_first(buffer, coordinates, 142298, 134918, 138158, 15,
                                             nmax);

        simdtrf::transform_f_inner(buffer, 144548, 142298, 15, 15, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 144548, 105, nmax);

        simdtrf::transform_f_inner(buffer, 144548, 140048, 15, 15, nmax);

        simdtrf::transform_g_outer(values + 945 * nvalues + n * npairs, nvalues, buffer, 144548,
                                   105, nmax);
    }

    for (size_t m = 0; m < 1890; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
