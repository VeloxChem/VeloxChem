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


#include "SimdThreeCenterElectronRepulsionRsRecDDK.hpp"

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
#include "SimdTransferDD.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_ddk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ddk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 34140, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 750 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 34140, 28008, 2937, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10, 11}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 18, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 16, 17,
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

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 81, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 84, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 87, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 90, 0, 3, 7, 8,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 96, 0, 3, 8, 9,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 9, 10,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 10, 11,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 114, 0, 3, 11, 12,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 120, 0, 3, 12, 13,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 126, 0, 3, 13, 14,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 14, 15,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 15, 16,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 144, 0, 3, 19, 20,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 150, 0, 3, 20, 21,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 156, 0, 3, 21, 22,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 22, 23,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 23, 24,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 174, 0, 3, 24, 25,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 180, 0, 3, 25, 26,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 186, 0, 3, 26, 27,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 27, 28,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 30, 33,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 33, 36,
                                                                       96, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 36, 39,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 39, 42,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 42, 45,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 45, 48,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 48, 51,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 51, 54,
                                                                       132, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 60, 63,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 63, 66,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 66, 69,
                                                                       156, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 69, 72,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 72, 75,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 75, 78,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 78, 81,
                                                                       180, 186, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 81, 84,
                                                                       186, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 90, 96,
                                                                       198, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 373, 0, 3, 96,
                                                                       102, 208, 218, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 102,
                                                                       108, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 403, 0, 3, 108,
                                                                       114, 228, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 114,
                                                                       120, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 433, 0, 3, 120,
                                                                       126, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 126,
                                                                       132, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 463, 0, 3, 144,
                                                                       150, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 150,
                                                                       156, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 493, 0, 3, 156,
                                                                       162, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 162,
                                                                       168, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 523, 0, 3, 168,
                                                                       174, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 538, 0, 3, 174,
                                                                       180, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 553, 0, 3, 180,
                                                                       186, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 568, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 571, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 574, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 577, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 580, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 583, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 586, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 589, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 592, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 595, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 598, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 601, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 604, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 607, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 610, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 613, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 616, 3, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 619, 3, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 622, 3, 26, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 625, 3, 27, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 628, 3, 28, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 631, 3, 29, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 634, 3, 7, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 643, 3, 8, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 652, 3, 9, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 661, 3, 10, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 670, 3, 11, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 679, 3, 12, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 688, 3, 13, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 697, 3, 14, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 706, 3, 15, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 715, 3, 16, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 724, 3, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 733, 3, 20, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 742, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 751, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 760, 3, 23, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 769, 3, 24, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 778, 3, 25, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 787, 3, 26, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 796, 3, 27, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 805, 3, 28, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 814, 3, 30, 90,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 832, 3, 33, 96,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 850, 3, 36, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 868, 3, 39, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 886, 3, 42, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 904, 3, 45, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 922, 3, 48, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 940, 3, 51, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 958, 3, 54, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 976, 3, 60, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 994, 3, 63, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1012, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1030, 3, 69, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1048, 3, 72, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1066, 3, 75, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1084, 3, 78, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1102, 3, 81, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1120, 3, 84, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1138, 3, 90, 198,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1168, 3, 96, 208,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1198, 3, 102, 218,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1228, 3, 108, 228,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1258, 3, 114, 238,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1288, 3, 120, 248,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1318, 3, 126, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1348, 3, 132, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1378, 3, 144, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1408, 3, 150, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1438, 3, 156, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1468, 3, 162, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1498, 3, 168, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1528, 3, 174, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1558, 3, 180, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1588, 3, 186, 348,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1618, 3, 198, 358,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1663, 3, 208, 373,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1708, 3, 218, 388,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1753, 3, 228, 403,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1798, 3, 238, 418,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1843, 3, 248, 433,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1888, 3, 258, 448,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1933, 3, 278, 463,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1978, 3, 288, 478,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2023, 3, 298, 493,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2068, 3, 308, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2113, 3, 318, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2158, 3, 328, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2203, 3, 338, 553,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2248, 3, 7, 8,
                                                                       574, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2254, 3, 8, 9,
                                                                       577, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2260, 3, 9, 10,
                                                                       580, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2266, 3, 10, 11,
                                                                       583, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2272, 3, 11, 12,
                                                                       586, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2278, 3, 12, 13,
                                                                       589, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2284, 3, 13, 14,
                                                                       592, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2290, 3, 14, 15,
                                                                       595, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2296, 3, 15, 16,
                                                                       598, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2302, 3, 19, 20,
                                                                       607, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2308, 3, 20, 21,
                                                                       610, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2314, 3, 21, 22,
                                                                       613, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2320, 3, 22, 23,
                                                                       616, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2326, 3, 23, 24,
                                                                       619, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2332, 3, 24, 25,
                                                                       622, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2338, 3, 25, 26,
                                                                       625, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2344, 3, 26, 27,
                                                                       628, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2350, 3, 27, 28,
                                                                       631, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2356, 0, 3, 2248,
                                                                       574, 2254, 652, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2374, 0, 3, 2254,
                                                                       577, 2260, 661, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2392, 0, 3, 2260,
                                                                       580, 2266, 670, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2410, 0, 3, 2266,
                                                                       583, 2272, 679, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2428, 0, 3, 2272,
                                                                       586, 2278, 688, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2446, 0, 3, 2278,
                                                                       589, 2284, 697, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2464, 0, 3, 2284,
                                                                       592, 2290, 706, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2482, 0, 3, 2290,
                                                                       595, 2296, 715, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2500, 0, 3, 2302,
                                                                       607, 2308, 742, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2518, 0, 3, 2308,
                                                                       610, 2314, 751, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2536, 0, 3, 2314,
                                                                       613, 2320, 760, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2554, 0, 3, 2320,
                                                                       616, 2326, 769, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2572, 0, 3, 2326,
                                                                       619, 2332, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2590, 0, 3, 2332,
                                                                       622, 2338, 787, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2608, 0, 3, 2338,
                                                                       625, 2344, 796, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2626, 0, 3, 2344,
                                                                       628, 2350, 805, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2644, 0, 3, 2356,
                                                                       652, 2374, 90, 96, 850,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2680, 0, 3, 2374,
                                                                       661, 2392, 96, 102, 868,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2716, 0, 3, 2392,
                                                                       670, 2410, 102, 108, 886,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2752, 0, 3, 2410,
                                                                       679, 2428, 108, 114, 904,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2788, 0, 3, 2428,
                                                                       688, 2446, 114, 120, 922,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2824, 0, 3, 2446,
                                                                       697, 2464, 120, 126, 940,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2860, 0, 3, 2464,
                                                                       706, 2482, 126, 132, 958,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2896, 0, 3, 2500,
                                                                       742, 2518, 144, 150, 1012,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2932, 0, 3, 2518,
                                                                       751, 2536, 150, 156, 1030,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2968, 0, 3, 2536,
                                                                       760, 2554, 156, 162, 1048,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3004, 0, 3, 2554,
                                                                       769, 2572, 162, 168, 1066,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3040, 0, 3, 2572,
                                                                       778, 2590, 168, 174, 1084,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3076, 0, 3, 2590,
                                                                       787, 2608, 174, 180, 1102,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3112, 0, 3, 2608,
                                                                       796, 2626, 180, 186, 1120,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3148, 0, 3, 2644,
                                                                       850, 2680, 198, 208, 1198,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3208, 0, 3, 2680,
                                                                       868, 2716, 208, 218, 1228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3268, 0, 3, 2716,
                                                                       886, 2752, 218, 228, 1258,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3328, 0, 3, 2752,
                                                                       904, 2788, 228, 238, 1288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3388, 0, 3, 2788,
                                                                       922, 2824, 238, 248, 1318,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3448, 0, 3, 2824,
                                                                       940, 2860, 248, 258, 1348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3508, 0, 3, 2896,
                                                                       1012, 2932, 278, 288,
                                                                       1438, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3568, 0, 3, 2932,
                                                                       1030, 2968, 288, 298,
                                                                       1468, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3628, 0, 3, 2968,
                                                                       1048, 3004, 298, 308,
                                                                       1498, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3688, 0, 3, 3004,
                                                                       1066, 3040, 308, 318,
                                                                       1528, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3748, 0, 3, 3040,
                                                                       1084, 3076, 318, 328,
                                                                       1558, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3808, 0, 3, 3076,
                                                                       1102, 3112, 328, 338,
                                                                       1588, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3868, 0, 3, 3148,
                                                                       1198, 3208, 358, 373,
                                                                       1708, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3958, 0, 3, 3208,
                                                                       1228, 3268, 373, 388,
                                                                       1753, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4048, 0, 3, 3268,
                                                                       1258, 3328, 388, 403,
                                                                       1798, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4138, 0, 3, 3328,
                                                                       1288, 3388, 403, 418,
                                                                       1843, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4228, 0, 3, 3388,
                                                                       1318, 3448, 418, 433,
                                                                       1888, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4318, 0, 3, 3508,
                                                                       1438, 3568, 463, 478,
                                                                       2023, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4408, 0, 3, 3568,
                                                                       1468, 3628, 478, 493,
                                                                       2068, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4498, 0, 3, 3628,
                                                                       1498, 3688, 493, 508,
                                                                       2113, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4588, 0, 3, 3688,
                                                                       1528, 3748, 508, 523,
                                                                       2158, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4678, 0, 3, 3748,
                                                                       1558, 3808, 523, 538,
                                                                       2203, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4768, 3, 568, 571,
                                                                       2248, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4778, 3, 571, 574,
                                                                       2254, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4788, 3, 574, 577,
                                                                       2260, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4798, 3, 577, 580,
                                                                       2266, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4808, 3, 580, 583,
                                                                       2272, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4818, 3, 583, 586,
                                                                       2278, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4828, 3, 586, 589,
                                                                       2284, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4838, 3, 589, 592,
                                                                       2290, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4848, 3, 592, 595,
                                                                       2296, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4858, 3, 601, 604,
                                                                       2302, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4868, 3, 604, 607,
                                                                       2308, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4878, 3, 607, 610,
                                                                       2314, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4888, 3, 610, 613,
                                                                       2320, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4898, 3, 613, 616,
                                                                       2326, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4908, 3, 616, 619,
                                                                       2332, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4918, 3, 619, 622,
                                                                       2338, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4928, 3, 622, 625,
                                                                       2344, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4938, 3, 625, 628,
                                                                       2350, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4948, 0, 3, 4768,
                                                                       2248, 4778, 634, 643,
                                                                       2356, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4978, 0, 3, 4778,
                                                                       2254, 4788, 643, 652,
                                                                       2374, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5008, 0, 3, 4788,
                                                                       2260, 4798, 652, 661,
                                                                       2392, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5038, 0, 3, 4798,
                                                                       2266, 4808, 661, 670,
                                                                       2410, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5068, 0, 3, 4808,
                                                                       2272, 4818, 670, 679,
                                                                       2428, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5098, 0, 3, 4818,
                                                                       2278, 4828, 679, 688,
                                                                       2446, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5128, 0, 3, 4828,
                                                                       2284, 4838, 688, 697,
                                                                       2464, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5158, 0, 3, 4838,
                                                                       2290, 4848, 697, 706,
                                                                       2482, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5188, 0, 3, 4858,
                                                                       2302, 4868, 724, 733,
                                                                       2500, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5218, 0, 3, 4868,
                                                                       2308, 4878, 733, 742,
                                                                       2518, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5248, 0, 3, 4878,
                                                                       2314, 4888, 742, 751,
                                                                       2536, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5278, 0, 3, 4888,
                                                                       2320, 4898, 751, 760,
                                                                       2554, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5308, 0, 3, 4898,
                                                                       2326, 4908, 760, 769,
                                                                       2572, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5338, 0, 3, 4908,
                                                                       2332, 4918, 769, 778,
                                                                       2590, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5368, 0, 3, 4918,
                                                                       2338, 4928, 778, 787,
                                                                       2608, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5398, 0, 3, 4928,
                                                                       2344, 4938, 787, 796,
                                                                       2626, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5428, 0, 3, 4948,
                                                                       2356, 4978, 814, 832,
                                                                       2644, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5488, 0, 3, 4978,
                                                                       2374, 5008, 832, 850,
                                                                       2680, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5548, 0, 3, 5008,
                                                                       2392, 5038, 850, 868,
                                                                       2716, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5608, 0, 3, 5038,
                                                                       2410, 5068, 868, 886,
                                                                       2752, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5668, 0, 3, 5068,
                                                                       2428, 5098, 886, 904,
                                                                       2788, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5728, 0, 3, 5098,
                                                                       2446, 5128, 904, 922,
                                                                       2824, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5788, 0, 3, 5128,
                                                                       2464, 5158, 922, 940,
                                                                       2860, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5848, 0, 3, 5188,
                                                                       2500, 5218, 976, 994,
                                                                       2896, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5908, 0, 3, 5218,
                                                                       2518, 5248, 994, 1012,
                                                                       2932, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5968, 0, 3, 5248,
                                                                       2536, 5278, 1012, 1030,
                                                                       2968, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6028, 0, 3, 5278,
                                                                       2554, 5308, 1030, 1048,
                                                                       3004, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6088, 0, 3, 5308,
                                                                       2572, 5338, 1048, 1066,
                                                                       3040, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6148, 0, 3, 5338,
                                                                       2590, 5368, 1066, 1084,
                                                                       3076, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6208, 0, 3, 5368,
                                                                       2608, 5398, 1084, 1102,
                                                                       3112, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6268, 0, 3, 5428,
                                                                       2644, 5488, 1138, 1168,
                                                                       3148, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6368, 0, 3, 5488,
                                                                       2680, 5548, 1168, 1198,
                                                                       3208, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6468, 0, 3, 5548,
                                                                       2716, 5608, 1198, 1228,
                                                                       3268, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6568, 0, 3, 5608,
                                                                       2752, 5668, 1228, 1258,
                                                                       3328, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6668, 0, 3, 5668,
                                                                       2788, 5728, 1258, 1288,
                                                                       3388, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6768, 0, 3, 5728,
                                                                       2824, 5788, 1288, 1318,
                                                                       3448, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6868, 0, 3, 5848,
                                                                       2896, 5908, 1378, 1408,
                                                                       3508, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6968, 0, 3, 5908,
                                                                       2932, 5968, 1408, 1438,
                                                                       3568, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7068, 0, 3, 5968,
                                                                       2968, 6028, 1438, 1468,
                                                                       3628, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7168, 0, 3, 6028,
                                                                       3004, 6088, 1468, 1498,
                                                                       3688, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7268, 0, 3, 6088,
                                                                       3040, 6148, 1498, 1528,
                                                                       3748, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7368, 0, 3, 6148,
                                                                       3076, 6208, 1528, 1558,
                                                                       3808, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7468, 0, 3, 6268,
                                                                       3148, 6368, 1618, 1663,
                                                                       3868, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7618, 0, 3, 6368,
                                                                       3208, 6468, 1663, 1708,
                                                                       3958, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7768, 0, 3, 6468,
                                                                       3268, 6568, 1708, 1753,
                                                                       4048, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7918, 0, 3, 6568,
                                                                       3328, 6668, 1753, 1798,
                                                                       4138, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8068, 0, 3, 6668,
                                                                       3388, 6768, 1798, 1843,
                                                                       4228, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8218, 0, 3, 6868,
                                                                       3508, 6968, 1933, 1978,
                                                                       4318, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8368, 0, 3, 6968,
                                                                       3568, 7068, 1978, 2023,
                                                                       4408, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8518, 0, 3, 7068,
                                                                       3628, 7168, 2023, 2068,
                                                                       4498, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8668, 0, 3, 7168,
                                                                       3688, 7268, 2068, 2113,
                                                                       4588, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8818, 0, 3, 7268,
                                                                       3748, 7368, 2113, 2158,
                                                                       4678, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8968, 3, 2248,
                                                                       2254, 4788, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8983, 3, 2254,
                                                                       2260, 4798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8998, 3, 2260,
                                                                       2266, 4808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9013, 3, 2266,
                                                                       2272, 4818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9028, 3, 2272,
                                                                       2278, 4828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9043, 3, 2278,
                                                                       2284, 4838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9058, 3, 2284,
                                                                       2290, 4848, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9073, 3, 2302,
                                                                       2308, 4878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9088, 3, 2308,
                                                                       2314, 4888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9103, 3, 2314,
                                                                       2320, 4898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9118, 3, 2320,
                                                                       2326, 4908, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9133, 3, 2326,
                                                                       2332, 4918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9148, 3, 2332,
                                                                       2338, 4928, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9163, 3, 2338,
                                                                       2344, 4938, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9178, 0, 3, 8968,
                                                                       4788, 8983, 2356, 2374,
                                                                       5008, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9223, 0, 3, 8983,
                                                                       4798, 8998, 2374, 2392,
                                                                       5038, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9268, 0, 3, 8998,
                                                                       4808, 9013, 2392, 2410,
                                                                       5068, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9313, 0, 3, 9013,
                                                                       4818, 9028, 2410, 2428,
                                                                       5098, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9358, 0, 3, 9028,
                                                                       4828, 9043, 2428, 2446,
                                                                       5128, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9403, 0, 3, 9043,
                                                                       4838, 9058, 2446, 2464,
                                                                       5158, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9448, 0, 3, 9073,
                                                                       4878, 9088, 2500, 2518,
                                                                       5248, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9493, 0, 3, 9088,
                                                                       4888, 9103, 2518, 2536,
                                                                       5278, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9538, 0, 3, 9103,
                                                                       4898, 9118, 2536, 2554,
                                                                       5308, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9583, 0, 3, 9118,
                                                                       4908, 9133, 2554, 2572,
                                                                       5338, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9628, 0, 3, 9133,
                                                                       4918, 9148, 2572, 2590,
                                                                       5368, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9673, 0, 3, 9148,
                                                                       4928, 9163, 2590, 2608,
                                                                       5398, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 9718, 0, 3, 9178,
                                                                       5008, 9223, 2644, 2680,
                                                                       5548, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 9808, 0, 3, 9223,
                                                                       5038, 9268, 2680, 2716,
                                                                       5608, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 9898, 0, 3, 9268,
                                                                       5068, 9313, 2716, 2752,
                                                                       5668, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 9988, 0, 3, 9313,
                                                                       5098, 9358, 2752, 2788,
                                                                       5728, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10078, 0, 3, 9358,
                                                                       5128, 9403, 2788, 2824,
                                                                       5788, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10168, 0, 3, 9448,
                                                                       5248, 9493, 2896, 2932,
                                                                       5968, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10258, 0, 3, 9493,
                                                                       5278, 9538, 2932, 2968,
                                                                       6028, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10348, 0, 3, 9538,
                                                                       5308, 9583, 2968, 3004,
                                                                       6088, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10438, 0, 3, 9583,
                                                                       5338, 9628, 3004, 3040,
                                                                       6148, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10528, 0, 3, 9628,
                                                                       5368, 9673, 3040, 3076,
                                                                       6208, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 10618, 0, 3, 9718,
                                                                       5548, 9808, 3148, 3208,
                                                                       6468, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 10768, 0, 3, 9808,
                                                                       5608, 9898, 3208, 3268,
                                                                       6568, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 10918, 0, 3, 9898,
                                                                       5668, 9988, 3268, 3328,
                                                                       6668, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11068, 0, 3, 9988,
                                                                       5728, 10078, 3328, 3388,
                                                                       6768, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11218, 0, 3,
                                                                       10168, 5968, 10258, 3508,
                                                                       3568, 7068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11368, 0, 3,
                                                                       10258, 6028, 10348, 3568,
                                                                       3628, 7168, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11518, 0, 3,
                                                                       10348, 6088, 10438, 3628,
                                                                       3688, 7268, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11668, 0, 3,
                                                                       10438, 6148, 10528, 3688,
                                                                       3748, 7368, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 11818, 0, 3,
                                                                       10618, 6468, 10768, 3868,
                                                                       3958, 7768, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 12043, 0, 3,
                                                                       10768, 6568, 10918, 3958,
                                                                       4048, 7918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 12268, 0, 3,
                                                                       10918, 6668, 11068, 4048,
                                                                       4138, 8068, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 12493, 0, 3,
                                                                       11218, 7068, 11368, 4318,
                                                                       4408, 8518, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 12718, 0, 3,
                                                                       11368, 7168, 11518, 4408,
                                                                       4498, 8668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 12943, 0, 3,
                                                                       11518, 7268, 11668, 4498,
                                                                       4588, 8818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13168, 3, 4768,
                                                                       4778, 8968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13189, 3, 4778,
                                                                       4788, 8983, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13210, 3, 4788,
                                                                       4798, 8998, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13231, 3, 4798,
                                                                       4808, 9013, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13252, 3, 4808,
                                                                       4818, 9028, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13273, 3, 4818,
                                                                       4828, 9043, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13294, 3, 4828,
                                                                       4838, 9058, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13315, 3, 4858,
                                                                       4868, 9073, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13336, 3, 4868,
                                                                       4878, 9088, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13357, 3, 4878,
                                                                       4888, 9103, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13378, 3, 4888,
                                                                       4898, 9118, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13399, 3, 4898,
                                                                       4908, 9133, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13420, 3, 4908,
                                                                       4918, 9148, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13441, 3, 4918,
                                                                       4928, 9163, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13462, 0, 3,
                                                                       13168, 8968, 13189, 4948,
                                                                       4978, 9178, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13525, 0, 3,
                                                                       13189, 8983, 13210, 4978,
                                                                       5008, 9223, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13588, 0, 3,
                                                                       13210, 8998, 13231, 5008,
                                                                       5038, 9268, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13651, 0, 3,
                                                                       13231, 9013, 13252, 5038,
                                                                       5068, 9313, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13714, 0, 3,
                                                                       13252, 9028, 13273, 5068,
                                                                       5098, 9358, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13777, 0, 3,
                                                                       13273, 9043, 13294, 5098,
                                                                       5128, 9403, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13840, 0, 3,
                                                                       13315, 9073, 13336, 5188,
                                                                       5218, 9448, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13903, 0, 3,
                                                                       13336, 9088, 13357, 5218,
                                                                       5248, 9493, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13966, 0, 3,
                                                                       13357, 9103, 13378, 5248,
                                                                       5278, 9538, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 14029, 0, 3,
                                                                       13378, 9118, 13399, 5278,
                                                                       5308, 9583, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 14092, 0, 3,
                                                                       13399, 9133, 13420, 5308,
                                                                       5338, 9628, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 14155, 0, 3,
                                                                       13420, 9148, 13441, 5338,
                                                                       5368, 9673, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14218, 0, 3,
                                                                       13462, 9178, 13525, 5428,
                                                                       5488, 9718, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14344, 0, 3,
                                                                       13525, 9223, 13588, 5488,
                                                                       5548, 9808, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14470, 0, 3,
                                                                       13588, 9268, 13651, 5548,
                                                                       5608, 9898, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14596, 0, 3,
                                                                       13651, 9313, 13714, 5608,
                                                                       5668, 9988, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14722, 0, 3,
                                                                       13714, 9358, 13777, 5668,
                                                                       5728, 10078, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14848, 0, 3,
                                                                       13840, 9448, 13903, 5848,
                                                                       5908, 10168, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14974, 0, 3,
                                                                       13903, 9493, 13966, 5908,
                                                                       5968, 10258, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 15100, 0, 3,
                                                                       13966, 9538, 14029, 5968,
                                                                       6028, 10348, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 15226, 0, 3,
                                                                       14029, 9583, 14092, 6028,
                                                                       6088, 10438, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 15352, 0, 3,
                                                                       14092, 9628, 14155, 6088,
                                                                       6148, 10528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 15478, 0, 3,
                                                                       14218, 9718, 14344, 6268,
                                                                       6368, 10618, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 15688, 0, 3,
                                                                       14344, 9808, 14470, 6368,
                                                                       6468, 10768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 15898, 0, 3,
                                                                       14470, 9898, 14596, 6468,
                                                                       6568, 10918, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 16108, 0, 3,
                                                                       14596, 9988, 14722, 6568,
                                                                       6668, 11068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 16318, 0, 3,
                                                                       14848, 10168, 14974, 6868,
                                                                       6968, 11218, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 16528, 0, 3,
                                                                       14974, 10258, 15100, 6968,
                                                                       7068, 11368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 16738, 0, 3,
                                                                       15100, 10348, 15226, 7068,
                                                                       7168, 11518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 16948, 0, 3,
                                                                       15226, 10438, 15352, 7168,
                                                                       7268, 11668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 17158, 0, 3,
                                                                       15478, 10618, 15688, 7468,
                                                                       7618, 11818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 17473, 0, 3,
                                                                       15688, 10768, 15898, 7618,
                                                                       7768, 12043, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 17788, 0, 3,
                                                                       15898, 10918, 16108, 7768,
                                                                       7918, 12268, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 18103, 0, 3,
                                                                       16318, 11218, 16528, 8218,
                                                                       8368, 12493, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 18418, 0, 3,
                                                                       16528, 11368, 16738, 8368,
                                                                       8518, 12718, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 18733, 0, 3,
                                                                       16738, 11518, 16948, 8518,
                                                                       8668, 12943, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19048, 3, 8968,
                                                                       8983, 13210, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19076, 3, 8983,
                                                                       8998, 13231, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19104, 3, 8998,
                                                                       9013, 13252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19132, 3, 9013,
                                                                       9028, 13273, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19160, 3, 9028,
                                                                       9043, 13294, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19188, 3, 9073,
                                                                       9088, 13357, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19216, 3, 9088,
                                                                       9103, 13378, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19244, 3, 9103,
                                                                       9118, 13399, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19272, 3, 9118,
                                                                       9133, 13420, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 19300, 3, 9133,
                                                                       9148, 13441, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19328, 0, 3,
                                                                       19048, 13210, 19076, 9178,
                                                                       9223, 13588, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19412, 0, 3,
                                                                       19076, 13231, 19104, 9223,
                                                                       9268, 13651, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19496, 0, 3,
                                                                       19104, 13252, 19132, 9268,
                                                                       9313, 13714, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19580, 0, 3,
                                                                       19132, 13273, 19160, 9313,
                                                                       9358, 13777, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19664, 0, 3,
                                                                       19188, 13357, 19216, 9448,
                                                                       9493, 13966, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19748, 0, 3,
                                                                       19216, 13378, 19244, 9493,
                                                                       9538, 14029, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19832, 0, 3,
                                                                       19244, 13399, 19272, 9538,
                                                                       9583, 14092, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 19916, 0, 3,
                                                                       19272, 13420, 19300, 9583,
                                                                       9628, 14155, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 20000, 0, 3,
                                                                       19328, 13588, 19412, 9718,
                                                                       9808, 14470, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 20168, 0, 3,
                                                                       19412, 13651, 19496, 9808,
                                                                       9898, 14596, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 20336, 0, 3,
                                                                       19496, 13714, 19580, 9898,
                                                                       9988, 14722, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 20504, 0, 3,
                                                                       19664, 13966, 19748,
                                                                       10168, 10258, 15100,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 20672, 0, 3,
                                                                       19748, 14029, 19832,
                                                                       10258, 10348, 15226,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 20840, 0, 3,
                                                                       19832, 14092, 19916,
                                                                       10348, 10438, 15352,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 21008, 0, 3,
                                                                       20000, 14470, 20168,
                                                                       10618, 10768, 15898,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 21288, 0, 3,
                                                                       20168, 14596, 20336,
                                                                       10768, 10918, 16108,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 21568, 0, 3,
                                                                       20504, 15100, 20672,
                                                                       11218, 11368, 16738,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 21848, 0, 3,
                                                                       20672, 15226, 20840,
                                                                       11368, 11518, 16948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 22128, 0, 3,
                                                                       21008, 15898, 21288,
                                                                       11818, 12043, 17788,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 22548, 0, 3,
                                                                       21568, 16738, 21848,
                                                                       12493, 12718, 18733,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 22968, 3, 13168,
                                                                       13189, 19048, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 23004, 3, 13189,
                                                                       13210, 19076, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 23040, 3, 13210,
                                                                       13231, 19104, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 23076, 3, 13231,
                                                                       13252, 19132, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 23112, 3, 13252,
                                                                       13273, 19160, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 23148, 3, 13315,
                                                                       13336, 19188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 23184, 3, 13336,
                                                                       13357, 19216, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 23220, 3, 13357,
                                                                       13378, 19244, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 23256, 3, 13378,
                                                                       13399, 19272, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 23292, 3, 13399,
                                                                       13420, 19300, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 23328, 0, 3,
                                                                       22968, 19048, 23004,
                                                                       13462, 13525, 19328,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 23436, 0, 3,
                                                                       23004, 19076, 23040,
                                                                       13525, 13588, 19412,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 23544, 0, 3,
                                                                       23040, 19104, 23076,
                                                                       13588, 13651, 19496,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 23652, 0, 3,
                                                                       23076, 19132, 23112,
                                                                       13651, 13714, 19580,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 23760, 0, 3,
                                                                       23148, 19188, 23184,
                                                                       13840, 13903, 19664,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 23868, 0, 3,
                                                                       23184, 19216, 23220,
                                                                       13903, 13966, 19748,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 23976, 0, 3,
                                                                       23220, 19244, 23256,
                                                                       13966, 14029, 19832,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 24084, 0, 3,
                                                                       23256, 19272, 23292,
                                                                       14029, 14092, 19916,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 24192, 0, 3,
                                                                       23328, 19328, 23436,
                                                                       14218, 14344, 20000,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 24408, 0, 3,
                                                                       23436, 19412, 23544,
                                                                       14344, 14470, 20168,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 24624, 0, 3,
                                                                       23544, 19496, 23652,
                                                                       14470, 14596, 20336,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 24840, 0, 3,
                                                                       23760, 19664, 23868,
                                                                       14848, 14974, 20504,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 25056, 0, 3,
                                                                       23868, 19748, 23976,
                                                                       14974, 15100, 20672,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 25272, 0, 3,
                                                                       23976, 19832, 24084,
                                                                       15100, 15226, 20840,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 25488, 0, 3,
                                                                       24192, 20000, 24408,
                                                                       15478, 15688, 21008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 25848, 0, 3,
                                                                       24408, 20168, 24624,
                                                                       15688, 15898, 21288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 26208, 0, 3,
                                                                       24840, 20504, 25056,
                                                                       16318, 16528, 21568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 26568, 0, 3,
                                                                       25056, 20672, 25272,
                                                                       16528, 16738, 21848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 26928, 0, 3,
                                                                       25488, 21008, 25848,
                                                                       17158, 17473, 22128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 27468, 0, 3,
                                                                       26208, 21568, 26568,
                                                                       18103, 18418, 22548,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 28008, 24192, 216, ncols);

                    simdfunc::contract_primitives(buffer, 28314, 24840, 216, ncols);

                    simdfunc::contract_primitives(buffer, 28620, 25488, 360, ncols);

                    simdfunc::contract_primitives(buffer, 29130, 26208, 360, ncols);

                    simdfunc::contract_primitives(buffer, 29640, 26928, 540, ncols);

                    simdfunc::contract_primitives(buffer, 30405, 27468, 540, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 28224, 28008, 6, 1, nmax);

        simdtrf::transform_k_inner(buffer, 28530, 28314, 6, 1, nmax);

        simdtrf::transform_k_inner(buffer, 28980, 28620, 10, 1, nmax);

        simdtrf::transform_k_inner(buffer, 29490, 29130, 10, 1, nmax);

        simdtrf::transform_k_inner(buffer, 30180, 29640, 15, 1, nmax);

        simdtrf::transform_k_inner(buffer, 30945, 30405, 15, 1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 31170, 28224, 28980, 15,
                                             nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 31440, 28530, 29490, 15,
                                             nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 31710, 28980, 30180, 15,
                                             nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 32160, 29490, 30945, 15,
                                             nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 32610, 31170, 31710, 15,
                                             nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 33150, 31440, 32160, 15,
                                             nmax);

        simdtrf::transform_d_inner(buffer, 33690, 33150, 6, 15, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 33690, 75, nmax);

        simdtrf::transform_d_inner(buffer, 33690, 32610, 6, 15, nmax);

        simdtrf::transform_d_outer(values + 375 * nvalues + n * npairs, nvalues, buffer, 33690,
                                   75, nmax);
    }

    for (size_t m = 0; m < 750; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
