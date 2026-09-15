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


#include "SimdThreeCenterElectronRepulsionRecDDK.hpp"

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
compute_ddk_three_center_electron_repulsion(double               *values,
                                            const size_t          npairs,
                                            const size_t          natoms,
                                            const CBasisFunction &a_function,
                                            const CBasisFunction &b_function,
                                            const CBasisFunction &c_function,
                                            const CSimdMatrix    &coordinates,
                                            const CSimdMatrix    &c_coordinates,
                                            CSimdMatrix          &buffer,
                                            const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ddk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 17299, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 375 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 17299, 14008, 1356, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fa = -b_exps[j] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pa(buffer, coordinates, 0, nmax, fa);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                simdfunc::compute_pair_exponent(buffer, coordinates, 6, nmax, mu);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 7, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11}, ncols, fj, 6,
                                                        fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 19, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 8, 9,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 9, 10,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 10, 11,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 11, 12,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 12, 13,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 13, 14,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 14, 15,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 91, 0, 3, 15, 16,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 16, 17,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 19, 22,
                                                                       49, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 22, 25,
                                                                       55, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 123, 0, 3, 25, 28,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 28, 31,
                                                                       67, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 31, 34,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 34, 37,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 37, 40,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 40, 43,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 49, 55,
                                                                       103, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 55, 61,
                                                                       113, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 61, 67,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 67, 73,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 73, 79,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 79, 85,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 85, 91,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 288, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 291, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 294, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 297, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 300, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 303, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 306, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 309, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 312, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 315, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 318, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 321, 3, 8, 19,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 330, 3, 9, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 339, 3, 10, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 348, 3, 11, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 357, 3, 12, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 366, 3, 13, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 375, 3, 14, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 384, 3, 15, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 393, 3, 16, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 402, 3, 17, 46,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 411, 3, 19, 49,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 429, 3, 22, 55,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 447, 3, 25, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 465, 3, 28, 67,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 483, 3, 31, 73,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 501, 3, 34, 79,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 519, 3, 37, 85,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 537, 3, 40, 91,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 555, 3, 43, 97,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 573, 3, 49, 103,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 603, 3, 55, 113,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 633, 3, 61, 123,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 663, 3, 67, 133,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 693, 3, 73, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 723, 3, 79, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 753, 3, 85, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 783, 3, 91, 173,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 813, 3, 103, 183,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 858, 3, 113, 198,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 903, 3, 123, 213,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 948, 3, 133, 228,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 993, 3, 143, 243,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1038, 3, 153, 258,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1083, 3, 163, 273,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1128, 3, 8, 9,
                                                                       294, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1134, 3, 9, 10,
                                                                       297, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1140, 3, 10, 11,
                                                                       300, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1146, 3, 11, 12,
                                                                       303, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1152, 3, 12, 13,
                                                                       306, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1158, 3, 13, 14,
                                                                       309, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1164, 3, 14, 15,
                                                                       312, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1170, 3, 15, 16,
                                                                       315, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1176, 3, 16, 17,
                                                                       318, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1182, 0, 3, 1128,
                                                                       294, 1134, 339, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1200, 0, 3, 1134,
                                                                       297, 1140, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1218, 0, 3, 1140,
                                                                       300, 1146, 357, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1236, 0, 3, 1146,
                                                                       303, 1152, 366, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1254, 0, 3, 1152,
                                                                       306, 1158, 375, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1272, 0, 3, 1158,
                                                                       309, 1164, 384, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1290, 0, 3, 1164,
                                                                       312, 1170, 393, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1308, 0, 3, 1170,
                                                                       315, 1176, 402, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1326, 0, 3, 1182,
                                                                       339, 1200, 49, 55, 447,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1362, 0, 3, 1200,
                                                                       348, 1218, 55, 61, 465,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1398, 0, 3, 1218,
                                                                       357, 1236, 61, 67, 483,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1434, 0, 3, 1236,
                                                                       366, 1254, 67, 73, 501,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1470, 0, 3, 1254,
                                                                       375, 1272, 73, 79, 519,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1506, 0, 3, 1272,
                                                                       384, 1290, 79, 85, 537,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1542, 0, 3, 1290,
                                                                       393, 1308, 85, 91, 555,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1578, 0, 3, 1326,
                                                                       447, 1362, 103, 113, 633,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1638, 0, 3, 1362,
                                                                       465, 1398, 113, 123, 663,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1698, 0, 3, 1398,
                                                                       483, 1434, 123, 133, 693,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1758, 0, 3, 1434,
                                                                       501, 1470, 133, 143, 723,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1818, 0, 3, 1470,
                                                                       519, 1506, 143, 153, 753,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1878, 0, 3, 1506,
                                                                       537, 1542, 153, 163, 783,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 1938, 0, 3, 1578,
                                                                       633, 1638, 183, 198, 903,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2028, 0, 3, 1638,
                                                                       663, 1698, 198, 213, 948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2118, 0, 3, 1698,
                                                                       693, 1758, 213, 228, 993,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2208, 0, 3, 1758,
                                                                       723, 1818, 228, 243, 1038,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2298, 0, 3, 1818,
                                                                       753, 1878, 243, 258, 1083,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2388, 3, 288, 291,
                                                                       1128, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2398, 3, 291, 294,
                                                                       1134, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2408, 3, 294, 297,
                                                                       1140, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2418, 3, 297, 300,
                                                                       1146, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2428, 3, 300, 303,
                                                                       1152, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2438, 3, 303, 306,
                                                                       1158, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2448, 3, 306, 309,
                                                                       1164, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2458, 3, 309, 312,
                                                                       1170, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2468, 3, 312, 315,
                                                                       1176, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2478, 0, 3, 2388,
                                                                       1128, 2398, 321, 330,
                                                                       1182, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2508, 0, 3, 2398,
                                                                       1134, 2408, 330, 339,
                                                                       1200, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2538, 0, 3, 2408,
                                                                       1140, 2418, 339, 348,
                                                                       1218, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2568, 0, 3, 2418,
                                                                       1146, 2428, 348, 357,
                                                                       1236, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2598, 0, 3, 2428,
                                                                       1152, 2438, 357, 366,
                                                                       1254, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2628, 0, 3, 2438,
                                                                       1158, 2448, 366, 375,
                                                                       1272, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2658, 0, 3, 2448,
                                                                       1164, 2458, 375, 384,
                                                                       1290, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2688, 0, 3, 2458,
                                                                       1170, 2468, 384, 393,
                                                                       1308, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2718, 0, 3, 2478,
                                                                       1182, 2508, 411, 429,
                                                                       1326, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2778, 0, 3, 2508,
                                                                       1200, 2538, 429, 447,
                                                                       1362, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2838, 0, 3, 2538,
                                                                       1218, 2568, 447, 465,
                                                                       1398, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2898, 0, 3, 2568,
                                                                       1236, 2598, 465, 483,
                                                                       1434, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2958, 0, 3, 2598,
                                                                       1254, 2628, 483, 501,
                                                                       1470, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3018, 0, 3, 2628,
                                                                       1272, 2658, 501, 519,
                                                                       1506, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3078, 0, 3, 2658,
                                                                       1290, 2688, 519, 537,
                                                                       1542, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3138, 0, 3, 2718,
                                                                       1326, 2778, 573, 603,
                                                                       1578, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3238, 0, 3, 2778,
                                                                       1362, 2838, 603, 633,
                                                                       1638, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3338, 0, 3, 2838,
                                                                       1398, 2898, 633, 663,
                                                                       1698, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3438, 0, 3, 2898,
                                                                       1434, 2958, 663, 693,
                                                                       1758, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3538, 0, 3, 2958,
                                                                       1470, 3018, 693, 723,
                                                                       1818, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3638, 0, 3, 3018,
                                                                       1506, 3078, 723, 753,
                                                                       1878, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 3738, 0, 3, 3138,
                                                                       1578, 3238, 813, 858,
                                                                       1938, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 3888, 0, 3, 3238,
                                                                       1638, 3338, 858, 903,
                                                                       2028, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4038, 0, 3, 3338,
                                                                       1698, 3438, 903, 948,
                                                                       2118, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4188, 0, 3, 3438,
                                                                       1758, 3538, 948, 993,
                                                                       2208, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4338, 0, 3, 3538,
                                                                       1818, 3638, 993, 1038,
                                                                       2298, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4488, 3, 1128,
                                                                       1134, 2408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4503, 3, 1134,
                                                                       1140, 2418, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4518, 3, 1140,
                                                                       1146, 2428, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4533, 3, 1146,
                                                                       1152, 2438, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4548, 3, 1152,
                                                                       1158, 2448, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4563, 3, 1158,
                                                                       1164, 2458, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4578, 3, 1164,
                                                                       1170, 2468, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4593, 0, 3, 4488,
                                                                       2408, 4503, 1182, 1200,
                                                                       2538, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4638, 0, 3, 4503,
                                                                       2418, 4518, 1200, 1218,
                                                                       2568, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4683, 0, 3, 4518,
                                                                       2428, 4533, 1218, 1236,
                                                                       2598, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4728, 0, 3, 4533,
                                                                       2438, 4548, 1236, 1254,
                                                                       2628, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4773, 0, 3, 4548,
                                                                       2448, 4563, 1254, 1272,
                                                                       2658, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4818, 0, 3, 4563,
                                                                       2458, 4578, 1272, 1290,
                                                                       2688, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 4863, 0, 3, 4593,
                                                                       2538, 4638, 1326, 1362,
                                                                       2838, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 4953, 0, 3, 4638,
                                                                       2568, 4683, 1362, 1398,
                                                                       2898, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5043, 0, 3, 4683,
                                                                       2598, 4728, 1398, 1434,
                                                                       2958, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5133, 0, 3, 4728,
                                                                       2628, 4773, 1434, 1470,
                                                                       3018, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5223, 0, 3, 4773,
                                                                       2658, 4818, 1470, 1506,
                                                                       3078, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 5313, 0, 3, 4863,
                                                                       2838, 4953, 1578, 1638,
                                                                       3338, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 5463, 0, 3, 4953,
                                                                       2898, 5043, 1638, 1698,
                                                                       3438, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 5613, 0, 3, 5043,
                                                                       2958, 5133, 1698, 1758,
                                                                       3538, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 5763, 0, 3, 5133,
                                                                       3018, 5223, 1758, 1818,
                                                                       3638, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 5913, 0, 3, 5313,
                                                                       3338, 5463, 1938, 2028,
                                                                       4038, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 6138, 0, 3, 5463,
                                                                       3438, 5613, 2028, 2118,
                                                                       4188, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 6363, 0, 3, 5613,
                                                                       3538, 5763, 2118, 2208,
                                                                       4338, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6588, 3, 2388,
                                                                       2398, 4488, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6609, 3, 2398,
                                                                       2408, 4503, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6630, 3, 2408,
                                                                       2418, 4518, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6651, 3, 2418,
                                                                       2428, 4533, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6672, 3, 2428,
                                                                       2438, 4548, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6693, 3, 2438,
                                                                       2448, 4563, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6714, 3, 2448,
                                                                       2458, 4578, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 6735, 0, 3, 6588,
                                                                       4488, 6609, 2478, 2508,
                                                                       4593, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 6798, 0, 3, 6609,
                                                                       4503, 6630, 2508, 2538,
                                                                       4638, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 6861, 0, 3, 6630,
                                                                       4518, 6651, 2538, 2568,
                                                                       4683, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 6924, 0, 3, 6651,
                                                                       4533, 6672, 2568, 2598,
                                                                       4728, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 6987, 0, 3, 6672,
                                                                       4548, 6693, 2598, 2628,
                                                                       4773, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7050, 0, 3, 6693,
                                                                       4563, 6714, 2628, 2658,
                                                                       4818, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 7113, 0, 3, 6735,
                                                                       4593, 6798, 2718, 2778,
                                                                       4863, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 7239, 0, 3, 6798,
                                                                       4638, 6861, 2778, 2838,
                                                                       4953, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 7365, 0, 3, 6861,
                                                                       4683, 6924, 2838, 2898,
                                                                       5043, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 7491, 0, 3, 6924,
                                                                       4728, 6987, 2898, 2958,
                                                                       5133, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 7617, 0, 3, 6987,
                                                                       4773, 7050, 2958, 3018,
                                                                       5223, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 7743, 0, 3, 7113,
                                                                       4863, 7239, 3138, 3238,
                                                                       5313, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 7953, 0, 3, 7239,
                                                                       4953, 7365, 3238, 3338,
                                                                       5463, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 8163, 0, 3, 7365,
                                                                       5043, 7491, 3338, 3438,
                                                                       5613, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 8373, 0, 3, 7491,
                                                                       5133, 7617, 3438, 3538,
                                                                       5763, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 8583, 0, 3, 7743,
                                                                       5313, 7953, 3738, 3888,
                                                                       5913, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 8898, 0, 3, 7953,
                                                                       5463, 8163, 3888, 4038,
                                                                       6138, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 9213, 0, 3, 8163,
                                                                       5613, 8373, 4038, 4188,
                                                                       6363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 9528, 3, 4488,
                                                                       4503, 6630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 9556, 3, 4503,
                                                                       4518, 6651, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 9584, 3, 4518,
                                                                       4533, 6672, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 9612, 3, 4533,
                                                                       4548, 6693, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 9640, 3, 4548,
                                                                       4563, 6714, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 9668, 0, 3, 9528,
                                                                       6630, 9556, 4593, 4638,
                                                                       6861, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 9752, 0, 3, 9556,
                                                                       6651, 9584, 4638, 4683,
                                                                       6924, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 9836, 0, 3, 9584,
                                                                       6672, 9612, 4683, 4728,
                                                                       6987, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 9920, 0, 3, 9612,
                                                                       6693, 9640, 4728, 4773,
                                                                       7050, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 10004, 0, 3, 9668,
                                                                       6861, 9752, 4863, 4953,
                                                                       7365, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 10172, 0, 3, 9752,
                                                                       6924, 9836, 4953, 5043,
                                                                       7491, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 10340, 0, 3, 9836,
                                                                       6987, 9920, 5043, 5133,
                                                                       7617, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 10508, 0, 3,
                                                                       10004, 7365, 10172, 5313,
                                                                       5463, 8163, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 10788, 0, 3,
                                                                       10172, 7491, 10340, 5463,
                                                                       5613, 8373, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 11068, 0, 3,
                                                                       10508, 8163, 10788, 5913,
                                                                       6138, 9213, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 11488, 3, 6588,
                                                                       6609, 9528, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 11524, 3, 6609,
                                                                       6630, 9556, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 11560, 3, 6630,
                                                                       6651, 9584, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 11596, 3, 6651,
                                                                       6672, 9612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 11632, 3, 6672,
                                                                       6693, 9640, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 11668, 0, 3,
                                                                       11488, 9528, 11524, 6735,
                                                                       6798, 9668, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 11776, 0, 3,
                                                                       11524, 9556, 11560, 6798,
                                                                       6861, 9752, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 11884, 0, 3,
                                                                       11560, 9584, 11596, 6861,
                                                                       6924, 9836, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 11992, 0, 3,
                                                                       11596, 9612, 11632, 6924,
                                                                       6987, 9920, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 12100, 0, 3,
                                                                       11668, 9668, 11776, 7113,
                                                                       7239, 10004, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 12316, 0, 3,
                                                                       11776, 9752, 11884, 7239,
                                                                       7365, 10172, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 12532, 0, 3,
                                                                       11884, 9836, 11992, 7365,
                                                                       7491, 10340, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 12748, 0, 3,
                                                                       12100, 10004, 12316, 7743,
                                                                       7953, 10508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 13108, 0, 3,
                                                                       12316, 10172, 12532, 7953,
                                                                       8163, 10788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 13468, 0, 3,
                                                                       12748, 10508, 13108, 8583,
                                                                       8898, 11068, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 14008, 12100, 216, ncols);

                    simdfunc::contract_primitives(buffer, 14314, 12748, 360, ncols);

                    simdfunc::contract_primitives(buffer, 14824, 13468, 540, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 14224, 14008, 6, 1, nmax);

        simdtrf::transform_k_inner(buffer, 14674, 14314, 10, 1, nmax);

        simdtrf::transform_k_inner(buffer, 15364, 14824, 15, 1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 15589, 14224, 14674, 15,
                                             nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 15859, 14674, 15364, 15,
                                             nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 16309, 15589, 15859, 15,
                                             nmax);

        simdtrf::transform_d_inner(buffer, 16849, 16309, 6, 15, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 16849, 75, nmax);
    }

    for (size_t m = 0; m < 375; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
