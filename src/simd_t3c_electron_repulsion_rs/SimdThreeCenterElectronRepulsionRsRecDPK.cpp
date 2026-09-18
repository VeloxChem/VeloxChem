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


#include "SimdThreeCenterElectronRepulsionRsRecDPK.hpp"

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
#include "SimdTransferDP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_dpk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_dpk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 17250, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 450 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 17250, 14808, 1482, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 17, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 7, 8,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 8, 9,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 9, 10,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 10, 11,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 11, 12,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 12, 13,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 13, 14,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 124, 0, 3, 14, 15,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 130, 0, 3, 18, 19,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 136, 0, 3, 19, 20,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 20, 21,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 21, 22,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 154, 0, 3, 22, 23,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 160, 0, 3, 23, 24,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 166, 0, 3, 24, 25,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 25, 26,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 28, 31,
                                                                       82, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 31, 34,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 34, 37,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 37, 40,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 40, 43,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 43, 46,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 46, 49,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 55, 58,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 58, 61,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 61, 64,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 64, 67,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 67, 70,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 70, 73,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 73, 76,
                                                                       166, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 318, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 321, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 324, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 327, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 330, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 333, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 336, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 339, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 342, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 345, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 348, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 351, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 354, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 357, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 360, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 363, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 366, 3, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 369, 3, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 372, 3, 26, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 375, 3, 27, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 378, 3, 7, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 387, 3, 8, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 396, 3, 9, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 405, 3, 10, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 414, 3, 11, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 423, 3, 12, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 432, 3, 13, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 441, 3, 14, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 450, 3, 15, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 459, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 468, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 477, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 486, 3, 21, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 495, 3, 22, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 504, 3, 23, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 513, 3, 24, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 522, 3, 25, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 531, 3, 26, 79,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 540, 3, 28, 82,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 558, 3, 31, 88,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 576, 3, 34, 94,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 594, 3, 37, 100,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 612, 3, 40, 106,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 630, 3, 43, 112,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 648, 3, 46, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 666, 3, 49, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 684, 3, 55, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 702, 3, 58, 136,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 720, 3, 61, 142,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 738, 3, 64, 148,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 756, 3, 67, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 774, 3, 70, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 792, 3, 73, 166,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 810, 3, 76, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 828, 3, 82, 178,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 858, 3, 88, 188,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 888, 3, 94, 198,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 918, 3, 100, 208,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 948, 3, 106, 218,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 978, 3, 112, 228,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1008, 3, 118, 238,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1038, 3, 130, 248,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1068, 3, 136, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1098, 3, 142, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1128, 3, 148, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1158, 3, 154, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1188, 3, 160, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1218, 3, 166, 308,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1248, 3, 7, 8,
                                                                       324, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1254, 3, 8, 9,
                                                                       327, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1260, 3, 9, 10,
                                                                       330, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1266, 3, 10, 11,
                                                                       333, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1272, 3, 11, 12,
                                                                       336, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1278, 3, 12, 13,
                                                                       339, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1284, 3, 13, 14,
                                                                       342, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1290, 3, 14, 15,
                                                                       345, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1296, 3, 18, 19,
                                                                       354, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1302, 3, 19, 20,
                                                                       357, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1308, 3, 20, 21,
                                                                       360, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1314, 3, 21, 22,
                                                                       363, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1320, 3, 22, 23,
                                                                       366, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1326, 3, 23, 24,
                                                                       369, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1332, 3, 24, 25,
                                                                       372, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1338, 3, 25, 26,
                                                                       375, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1344, 0, 3, 1248,
                                                                       324, 1254, 396, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1362, 0, 3, 1254,
                                                                       327, 1260, 405, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 1260,
                                                                       330, 1266, 414, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1398, 0, 3, 1266,
                                                                       333, 1272, 423, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1416, 0, 3, 1272,
                                                                       336, 1278, 432, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1434, 0, 3, 1278,
                                                                       339, 1284, 441, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1452, 0, 3, 1284,
                                                                       342, 1290, 450, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1470, 0, 3, 1296,
                                                                       354, 1302, 477, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1488, 0, 3, 1302,
                                                                       357, 1308, 486, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1506, 0, 3, 1308,
                                                                       360, 1314, 495, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1524, 0, 3, 1314,
                                                                       363, 1320, 504, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1542, 0, 3, 1320,
                                                                       366, 1326, 513, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1560, 0, 3, 1326,
                                                                       369, 1332, 522, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1578, 0, 3, 1332,
                                                                       372, 1338, 531, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1596, 0, 3, 1344,
                                                                       396, 1362, 82, 88, 576,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 1362,
                                                                       405, 1380, 88, 94, 594,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1668, 0, 3, 1380,
                                                                       414, 1398, 94, 100, 612,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1704, 0, 3, 1398,
                                                                       423, 1416, 100, 106, 630,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1740, 0, 3, 1416,
                                                                       432, 1434, 106, 112, 648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1776, 0, 3, 1434,
                                                                       441, 1452, 112, 118, 666,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1812, 0, 3, 1470,
                                                                       477, 1488, 130, 136, 720,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1848, 0, 3, 1488,
                                                                       486, 1506, 136, 142, 738,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 1506,
                                                                       495, 1524, 142, 148, 756,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1920, 0, 3, 1524,
                                                                       504, 1542, 148, 154, 774,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1956, 0, 3, 1542,
                                                                       513, 1560, 154, 160, 792,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1992, 0, 3, 1560,
                                                                       522, 1578, 160, 166, 810,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2028, 0, 3, 1596,
                                                                       576, 1632, 178, 188, 888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2088, 0, 3, 1632,
                                                                       594, 1668, 188, 198, 918,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2148, 0, 3, 1668,
                                                                       612, 1704, 198, 208, 948,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2208, 0, 3, 1704,
                                                                       630, 1740, 208, 218, 978,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2268, 0, 3, 1740,
                                                                       648, 1776, 218, 228, 1008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2328, 0, 3, 1812,
                                                                       720, 1848, 248, 258, 1098,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2388, 0, 3, 1848,
                                                                       738, 1884, 258, 268, 1128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2448, 0, 3, 1884,
                                                                       756, 1920, 268, 278, 1158,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2508, 0, 3, 1920,
                                                                       774, 1956, 278, 288, 1188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2568, 0, 3, 1956,
                                                                       792, 1992, 288, 298, 1218,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2628, 3, 318, 321,
                                                                       1248, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2638, 3, 321, 324,
                                                                       1254, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2648, 3, 324, 327,
                                                                       1260, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2658, 3, 327, 330,
                                                                       1266, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2668, 3, 330, 333,
                                                                       1272, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2678, 3, 333, 336,
                                                                       1278, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2688, 3, 336, 339,
                                                                       1284, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2698, 3, 339, 342,
                                                                       1290, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2708, 3, 348, 351,
                                                                       1296, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2718, 3, 351, 354,
                                                                       1302, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2728, 3, 354, 357,
                                                                       1308, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2738, 3, 357, 360,
                                                                       1314, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2748, 3, 360, 363,
                                                                       1320, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2758, 3, 363, 366,
                                                                       1326, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2768, 3, 366, 369,
                                                                       1332, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2778, 3, 369, 372,
                                                                       1338, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2788, 0, 3, 2628,
                                                                       1248, 2638, 378, 387,
                                                                       1344, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2818, 0, 3, 2638,
                                                                       1254, 2648, 387, 396,
                                                                       1362, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2848, 0, 3, 2648,
                                                                       1260, 2658, 396, 405,
                                                                       1380, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2878, 0, 3, 2658,
                                                                       1266, 2668, 405, 414,
                                                                       1398, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2908, 0, 3, 2668,
                                                                       1272, 2678, 414, 423,
                                                                       1416, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2938, 0, 3, 2678,
                                                                       1278, 2688, 423, 432,
                                                                       1434, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2968, 0, 3, 2688,
                                                                       1284, 2698, 432, 441,
                                                                       1452, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2998, 0, 3, 2708,
                                                                       1296, 2718, 459, 468,
                                                                       1470, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3028, 0, 3, 2718,
                                                                       1302, 2728, 468, 477,
                                                                       1488, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3058, 0, 3, 2728,
                                                                       1308, 2738, 477, 486,
                                                                       1506, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 2738,
                                                                       1314, 2748, 486, 495,
                                                                       1524, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3118, 0, 3, 2748,
                                                                       1320, 2758, 495, 504,
                                                                       1542, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3148, 0, 3, 2758,
                                                                       1326, 2768, 504, 513,
                                                                       1560, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3178, 0, 3, 2768,
                                                                       1332, 2778, 513, 522,
                                                                       1578, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3208, 0, 3, 2788,
                                                                       1344, 2818, 540, 558,
                                                                       1596, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3268, 0, 3, 2818,
                                                                       1362, 2848, 558, 576,
                                                                       1632, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3328, 0, 3, 2848,
                                                                       1380, 2878, 576, 594,
                                                                       1668, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3388, 0, 3, 2878,
                                                                       1398, 2908, 594, 612,
                                                                       1704, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3448, 0, 3, 2908,
                                                                       1416, 2938, 612, 630,
                                                                       1740, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3508, 0, 3, 2938,
                                                                       1434, 2968, 630, 648,
                                                                       1776, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3568, 0, 3, 2998,
                                                                       1470, 3028, 684, 702,
                                                                       1812, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3628, 0, 3, 3028,
                                                                       1488, 3058, 702, 720,
                                                                       1848, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3688, 0, 3, 3058,
                                                                       1506, 3088, 720, 738,
                                                                       1884, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3748, 0, 3, 3088,
                                                                       1524, 3118, 738, 756,
                                                                       1920, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3808, 0, 3, 3118,
                                                                       1542, 3148, 756, 774,
                                                                       1956, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3868, 0, 3, 3148,
                                                                       1560, 3178, 774, 792,
                                                                       1992, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3928, 0, 3, 3208,
                                                                       1596, 3268, 828, 858,
                                                                       2028, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4028, 0, 3, 3268,
                                                                       1632, 3328, 858, 888,
                                                                       2088, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4128, 0, 3, 3328,
                                                                       1668, 3388, 888, 918,
                                                                       2148, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4228, 0, 3, 3388,
                                                                       1704, 3448, 918, 948,
                                                                       2208, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4328, 0, 3, 3448,
                                                                       1740, 3508, 948, 978,
                                                                       2268, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4428, 0, 3, 3568,
                                                                       1812, 3628, 1038, 1068,
                                                                       2328, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4528, 0, 3, 3628,
                                                                       1848, 3688, 1068, 1098,
                                                                       2388, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4628, 0, 3, 3688,
                                                                       1884, 3748, 1098, 1128,
                                                                       2448, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4728, 0, 3, 3748,
                                                                       1920, 3808, 1128, 1158,
                                                                       2508, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4828, 0, 3, 3808,
                                                                       1956, 3868, 1158, 1188,
                                                                       2568, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4928, 3, 1248,
                                                                       1254, 2648, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4943, 3, 1254,
                                                                       1260, 2658, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4958, 3, 1260,
                                                                       1266, 2668, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4973, 3, 1266,
                                                                       1272, 2678, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4988, 3, 1272,
                                                                       1278, 2688, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5003, 3, 1278,
                                                                       1284, 2698, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5018, 3, 1296,
                                                                       1302, 2728, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5033, 3, 1302,
                                                                       1308, 2738, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5048, 3, 1308,
                                                                       1314, 2748, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5063, 3, 1314,
                                                                       1320, 2758, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5078, 3, 1320,
                                                                       1326, 2768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5093, 3, 1326,
                                                                       1332, 2778, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5108, 0, 3, 4928,
                                                                       2648, 4943, 1344, 1362,
                                                                       2848, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5153, 0, 3, 4943,
                                                                       2658, 4958, 1362, 1380,
                                                                       2878, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5198, 0, 3, 4958,
                                                                       2668, 4973, 1380, 1398,
                                                                       2908, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5243, 0, 3, 4973,
                                                                       2678, 4988, 1398, 1416,
                                                                       2938, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 4988,
                                                                       2688, 5003, 1416, 1434,
                                                                       2968, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5333, 0, 3, 5018,
                                                                       2728, 5033, 1470, 1488,
                                                                       3058, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5378, 0, 3, 5033,
                                                                       2738, 5048, 1488, 1506,
                                                                       3088, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5423, 0, 3, 5048,
                                                                       2748, 5063, 1506, 1524,
                                                                       3118, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5468, 0, 3, 5063,
                                                                       2758, 5078, 1524, 1542,
                                                                       3148, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5513, 0, 3, 5078,
                                                                       2768, 5093, 1542, 1560,
                                                                       3178, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5558, 0, 3, 5108,
                                                                       2848, 5153, 1596, 1632,
                                                                       3328, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5648, 0, 3, 5153,
                                                                       2878, 5198, 1632, 1668,
                                                                       3388, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5738, 0, 3, 5198,
                                                                       2908, 5243, 1668, 1704,
                                                                       3448, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5828, 0, 3, 5243,
                                                                       2938, 5288, 1704, 1740,
                                                                       3508, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5918, 0, 3, 5333,
                                                                       3058, 5378, 1812, 1848,
                                                                       3688, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 6008, 0, 3, 5378,
                                                                       3088, 5423, 1848, 1884,
                                                                       3748, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 6098, 0, 3, 5423,
                                                                       3118, 5468, 1884, 1920,
                                                                       3808, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 6188, 0, 3, 5468,
                                                                       3148, 5513, 1920, 1956,
                                                                       3868, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6278, 0, 3, 5558,
                                                                       3328, 5648, 2028, 2088,
                                                                       4128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6428, 0, 3, 5648,
                                                                       3388, 5738, 2088, 2148,
                                                                       4228, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6578, 0, 3, 5738,
                                                                       3448, 5828, 2148, 2208,
                                                                       4328, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6728, 0, 3, 5918,
                                                                       3688, 6008, 2328, 2388,
                                                                       4628, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6878, 0, 3, 6008,
                                                                       3748, 6098, 2388, 2448,
                                                                       4728, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 7028, 0, 3, 6098,
                                                                       3808, 6188, 2448, 2508,
                                                                       4828, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7178, 3, 2628,
                                                                       2638, 4928, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7199, 3, 2638,
                                                                       2648, 4943, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7220, 3, 2648,
                                                                       2658, 4958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7241, 3, 2658,
                                                                       2668, 4973, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7262, 3, 2668,
                                                                       2678, 4988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7283, 3, 2678,
                                                                       2688, 5003, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7304, 3, 2708,
                                                                       2718, 5018, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7325, 3, 2718,
                                                                       2728, 5033, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7346, 3, 2728,
                                                                       2738, 5048, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7367, 3, 2738,
                                                                       2748, 5063, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7388, 3, 2748,
                                                                       2758, 5078, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7409, 3, 2758,
                                                                       2768, 5093, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7430, 0, 3, 7178,
                                                                       4928, 7199, 2788, 2818,
                                                                       5108, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7493, 0, 3, 7199,
                                                                       4943, 7220, 2818, 2848,
                                                                       5153, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7556, 0, 3, 7220,
                                                                       4958, 7241, 2848, 2878,
                                                                       5198, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7619, 0, 3, 7241,
                                                                       4973, 7262, 2878, 2908,
                                                                       5243, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7682, 0, 3, 7262,
                                                                       4988, 7283, 2908, 2938,
                                                                       5288, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7745, 0, 3, 7304,
                                                                       5018, 7325, 2998, 3028,
                                                                       5333, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7808, 0, 3, 7325,
                                                                       5033, 7346, 3028, 3058,
                                                                       5378, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7871, 0, 3, 7346,
                                                                       5048, 7367, 3058, 3088,
                                                                       5423, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7934, 0, 3, 7367,
                                                                       5063, 7388, 3088, 3118,
                                                                       5468, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7997, 0, 3, 7388,
                                                                       5078, 7409, 3118, 3148,
                                                                       5513, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 8060, 0, 3, 7430,
                                                                       5108, 7493, 3208, 3268,
                                                                       5558, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 8186, 0, 3, 7493,
                                                                       5153, 7556, 3268, 3328,
                                                                       5648, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 8312, 0, 3, 7556,
                                                                       5198, 7619, 3328, 3388,
                                                                       5738, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 8438, 0, 3, 7619,
                                                                       5243, 7682, 3388, 3448,
                                                                       5828, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 8564, 0, 3, 7745,
                                                                       5333, 7808, 3568, 3628,
                                                                       5918, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 8690, 0, 3, 7808,
                                                                       5378, 7871, 3628, 3688,
                                                                       6008, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 8816, 0, 3, 7871,
                                                                       5423, 7934, 3688, 3748,
                                                                       6098, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 8942, 0, 3, 7934,
                                                                       5468, 7997, 3748, 3808,
                                                                       6188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 9068, 0, 3, 8060,
                                                                       5558, 8186, 3928, 4028,
                                                                       6278, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 9278, 0, 3, 8186,
                                                                       5648, 8312, 4028, 4128,
                                                                       6428, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 9488, 0, 3, 8312,
                                                                       5738, 8438, 4128, 4228,
                                                                       6578, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 9698, 0, 3, 8564,
                                                                       5918, 8690, 4428, 4528,
                                                                       6728, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 9908, 0, 3, 8690,
                                                                       6008, 8816, 4528, 4628,
                                                                       6878, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 10118, 0, 3, 8816,
                                                                       6098, 8942, 4628, 4728,
                                                                       7028, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10328, 3, 4928,
                                                                       4943, 7220, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10356, 3, 4943,
                                                                       4958, 7241, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10384, 3, 4958,
                                                                       4973, 7262, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10412, 3, 4973,
                                                                       4988, 7283, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10440, 3, 5018,
                                                                       5033, 7346, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10468, 3, 5033,
                                                                       5048, 7367, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10496, 3, 5048,
                                                                       5063, 7388, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10524, 3, 5063,
                                                                       5078, 7409, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 10552, 0, 3,
                                                                       10328, 7220, 10356, 5108,
                                                                       5153, 7556, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 10636, 0, 3,
                                                                       10356, 7241, 10384, 5153,
                                                                       5198, 7619, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 10720, 0, 3,
                                                                       10384, 7262, 10412, 5198,
                                                                       5243, 7682, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 10804, 0, 3,
                                                                       10440, 7346, 10468, 5333,
                                                                       5378, 7871, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 10888, 0, 3,
                                                                       10468, 7367, 10496, 5378,
                                                                       5423, 7934, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 10972, 0, 3,
                                                                       10496, 7388, 10524, 5423,
                                                                       5468, 7997, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 11056, 0, 3,
                                                                       10552, 7556, 10636, 5558,
                                                                       5648, 8312, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 11224, 0, 3,
                                                                       10636, 7619, 10720, 5648,
                                                                       5738, 8438, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 11392, 0, 3,
                                                                       10804, 7871, 10888, 5918,
                                                                       6008, 8816, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 11560, 0, 3,
                                                                       10888, 7934, 10972, 6008,
                                                                       6098, 8942, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 11728, 0, 3,
                                                                       11056, 8312, 11224, 6278,
                                                                       6428, 9488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 12008, 0, 3,
                                                                       11392, 8816, 11560, 6728,
                                                                       6878, 10118, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12288, 3, 7178,
                                                                       7199, 10328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12324, 3, 7199,
                                                                       7220, 10356, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12360, 3, 7220,
                                                                       7241, 10384, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12396, 3, 7241,
                                                                       7262, 10412, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12432, 3, 7304,
                                                                       7325, 10440, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12468, 3, 7325,
                                                                       7346, 10468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12504, 3, 7346,
                                                                       7367, 10496, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12540, 3, 7367,
                                                                       7388, 10524, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 12576, 0, 3,
                                                                       12288, 10328, 12324, 7430,
                                                                       7493, 10552, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 12684, 0, 3,
                                                                       12324, 10356, 12360, 7493,
                                                                       7556, 10636, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 12792, 0, 3,
                                                                       12360, 10384, 12396, 7556,
                                                                       7619, 10720, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 12900, 0, 3,
                                                                       12432, 10440, 12468, 7745,
                                                                       7808, 10804, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 13008, 0, 3,
                                                                       12468, 10468, 12504, 7808,
                                                                       7871, 10888, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 13116, 0, 3,
                                                                       12504, 10496, 12540, 7871,
                                                                       7934, 10972, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 13224, 0, 3,
                                                                       12576, 10552, 12684, 8060,
                                                                       8186, 11056, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 13440, 0, 3,
                                                                       12684, 10636, 12792, 8186,
                                                                       8312, 11224, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 13656, 0, 3,
                                                                       12900, 10804, 13008, 8564,
                                                                       8690, 11392, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 13872, 0, 3,
                                                                       13008, 10888, 13116, 8690,
                                                                       8816, 11560, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 14088, 0, 3,
                                                                       13224, 11056, 13440, 9068,
                                                                       9278, 11728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 14448, 0, 3,
                                                                       13656, 11392, 13872, 9698,
                                                                       9908, 12008, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 14808, 13224, 216, ncols);

                    simdfunc::contract_primitives(buffer, 15114, 13656, 216, ncols);

                    simdfunc::contract_primitives(buffer, 15420, 14088, 360, ncols);

                    simdfunc::contract_primitives(buffer, 15930, 14448, 360, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 15024, 14808, 6, 1, nmax);

        simdtrf::transform_k_inner(buffer, 15330, 15114, 6, 1, nmax);

        simdtrf::transform_k_inner(buffer, 15780, 15420, 10, 1, nmax);

        simdtrf::transform_k_inner(buffer, 16290, 15930, 10, 1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 16440, 15024, 15780, 15,
                                             nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 16710, 15330, 16290, 15,
                                             nmax);

        simdtrf::transform_p_inner(buffer, 16980, 16710, 6, 15, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 16980, 45, nmax);

        simdtrf::transform_p_inner(buffer, 16980, 16440, 6, 15, nmax);

        simdtrf::transform_d_outer(values + 225 * nvalues + n * npairs, nvalues, buffer, 16980,
                                   45, nmax);
    }

    for (size_t m = 0; m < 450; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
