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


#include "SimdThreeCenterElectronRepulsionRecIDF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_idf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_idf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 16581, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 455 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 16581, 11228, 1538, dimensions);

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

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 103,
                                                                       113, 183, 198, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 309, 0, 3, 113,
                                                                       123, 198, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 330, 0, 3, 123,
                                                                       133, 213, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 351, 0, 3, 133,
                                                                       143, 228, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 372, 0, 3, 143,
                                                                       153, 243, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 393, 0, 3, 153,
                                                                       163, 258, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 414, 0, 3, 183,
                                                                       198, 288, 309, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 442, 0, 3, 198,
                                                                       213, 309, 330, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 470, 0, 3, 213,
                                                                       228, 330, 351, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 228,
                                                                       243, 351, 372, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 526, 0, 3, 243,
                                                                       258, 372, 393, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 554, 0, 3, 288,
                                                                       309, 414, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 590, 0, 3, 309,
                                                                       330, 442, 470, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 626, 0, 3, 330,
                                                                       351, 470, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 662, 0, 3, 351,
                                                                       372, 498, 526, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 698, 0, 3, 414,
                                                                       442, 554, 590, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 743, 0, 3, 442,
                                                                       470, 590, 626, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 788, 0, 3, 470,
                                                                       498, 626, 662, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 833, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 836, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 839, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 842, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 845, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 848, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 851, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 854, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 857, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 860, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 863, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 866, 3, 8, 19,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 875, 3, 9, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 884, 3, 10, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 893, 3, 11, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 902, 3, 12, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 911, 3, 13, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 920, 3, 14, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 929, 3, 15, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 938, 3, 16, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 947, 3, 17, 46,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 956, 3, 19, 49,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 974, 3, 22, 55,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 992, 3, 25, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1010, 3, 28, 67,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1028, 3, 31, 73,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1046, 3, 34, 79,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1064, 3, 37, 85,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1082, 3, 40, 91,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1100, 3, 43, 97,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1118, 3, 49, 103,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1148, 3, 55, 113,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1178, 3, 61, 123,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1208, 3, 67, 133,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1238, 3, 73, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1268, 3, 79, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1298, 3, 85, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1328, 3, 91, 173,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1358, 3, 103, 183,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1403, 3, 113, 198,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1448, 3, 123, 213,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1493, 3, 133, 228,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1538, 3, 143, 243,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1583, 3, 153, 258,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1628, 3, 163, 273,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1673, 3, 183, 288,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1736, 3, 198, 309,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1799, 3, 213, 330,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1862, 3, 228, 351,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1925, 3, 243, 372,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1988, 3, 258, 393,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2051, 3, 288, 414,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2135, 3, 309, 442,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2219, 3, 330, 470,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2303, 3, 351, 498,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2387, 3, 372, 526,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 2471, 3, 414, 554,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 2579, 3, 442, 590,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 2687, 3, 470, 626,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 2795, 3, 498, 662,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 2903, 3, 554, 698,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 3038, 3, 590, 743,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 3173, 3, 626, 788,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3308, 3, 8, 9,
                                                                       839, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3314, 3, 9, 10,
                                                                       842, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3320, 3, 10, 11,
                                                                       845, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3326, 3, 11, 12,
                                                                       848, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3332, 3, 12, 13,
                                                                       851, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3338, 3, 13, 14,
                                                                       854, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3344, 3, 14, 15,
                                                                       857, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3350, 3, 15, 16,
                                                                       860, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3356, 3, 16, 17,
                                                                       863, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3362, 0, 3, 3308,
                                                                       839, 3314, 884, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3380, 0, 3, 3314,
                                                                       842, 3320, 893, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3398, 0, 3, 3320,
                                                                       845, 3326, 902, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3416, 0, 3, 3326,
                                                                       848, 3332, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3434, 0, 3, 3332,
                                                                       851, 3338, 920, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3452, 0, 3, 3338,
                                                                       854, 3344, 929, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3470, 0, 3, 3344,
                                                                       857, 3350, 938, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3488, 0, 3, 3350,
                                                                       860, 3356, 947, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3506, 0, 3, 3362,
                                                                       884, 3380, 49, 55, 992,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3542, 0, 3, 3380,
                                                                       893, 3398, 55, 61, 1010,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3578, 0, 3, 3398,
                                                                       902, 3416, 61, 67, 1028,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3614, 0, 3, 3416,
                                                                       911, 3434, 67, 73, 1046,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3650, 0, 3, 3434,
                                                                       920, 3452, 73, 79, 1064,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3686, 0, 3, 3452,
                                                                       929, 3470, 79, 85, 1082,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3722, 0, 3, 3470,
                                                                       938, 3488, 85, 91, 1100,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3758, 0, 3, 3506,
                                                                       992, 3542, 103, 113, 1178,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3818, 0, 3, 3542,
                                                                       1010, 3578, 113, 123,
                                                                       1208, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3878, 0, 3, 3578,
                                                                       1028, 3614, 123, 133,
                                                                       1238, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3938, 0, 3, 3614,
                                                                       1046, 3650, 133, 143,
                                                                       1268, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3998, 0, 3, 3650,
                                                                       1064, 3686, 143, 153,
                                                                       1298, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4058, 0, 3, 3686,
                                                                       1082, 3722, 153, 163,
                                                                       1328, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4118, 0, 3, 3758,
                                                                       1178, 3818, 183, 198,
                                                                       1448, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4208, 0, 3, 3818,
                                                                       1208, 3878, 198, 213,
                                                                       1493, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 3878,
                                                                       1238, 3938, 213, 228,
                                                                       1538, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4388, 0, 3, 3938,
                                                                       1268, 3998, 228, 243,
                                                                       1583, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4478, 0, 3, 3998,
                                                                       1298, 4058, 243, 258,
                                                                       1628, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4568, 0, 3, 4118,
                                                                       1448, 4208, 288, 309,
                                                                       1799, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4694, 0, 3, 4208,
                                                                       1493, 4298, 309, 330,
                                                                       1862, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4820, 0, 3, 4298,
                                                                       1538, 4388, 330, 351,
                                                                       1925, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4946, 0, 3, 4388,
                                                                       1583, 4478, 351, 372,
                                                                       1988, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5072, 0, 3, 4568,
                                                                       1799, 4694, 414, 442,
                                                                       2219, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5240, 0, 3, 4694,
                                                                       1862, 4820, 442, 470,
                                                                       2303, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5408, 0, 3, 4820,
                                                                       1925, 4946, 470, 498,
                                                                       2387, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 5576, 0, 3, 5072,
                                                                       2219, 5240, 554, 590,
                                                                       2687, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 5792, 0, 3, 5240,
                                                                       2303, 5408, 590, 626,
                                                                       2795, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 6008, 0, 3, 5576,
                                                                       2687, 5792, 698, 743,
                                                                       3173, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6278, 3, 833, 836,
                                                                       3308, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6288, 3, 836, 839,
                                                                       3314, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6298, 3, 839, 842,
                                                                       3320, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6308, 3, 842, 845,
                                                                       3326, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6318, 3, 845, 848,
                                                                       3332, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6328, 3, 848, 851,
                                                                       3338, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6338, 3, 851, 854,
                                                                       3344, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6348, 3, 854, 857,
                                                                       3350, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6358, 3, 857, 860,
                                                                       3356, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6368, 0, 3, 6278,
                                                                       3308, 6288, 866, 875,
                                                                       3362, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6398, 0, 3, 6288,
                                                                       3314, 6298, 875, 884,
                                                                       3380, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6428, 0, 3, 6298,
                                                                       3320, 6308, 884, 893,
                                                                       3398, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6458, 0, 3, 6308,
                                                                       3326, 6318, 893, 902,
                                                                       3416, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6488, 0, 3, 6318,
                                                                       3332, 6328, 902, 911,
                                                                       3434, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6518, 0, 3, 6328,
                                                                       3338, 6338, 911, 920,
                                                                       3452, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6548, 0, 3, 6338,
                                                                       3344, 6348, 920, 929,
                                                                       3470, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6578, 0, 3, 6348,
                                                                       3350, 6358, 929, 938,
                                                                       3488, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 6368,
                                                                       3362, 6398, 956, 974,
                                                                       3506, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6668, 0, 3, 6398,
                                                                       3380, 6428, 974, 992,
                                                                       3542, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6728, 0, 3, 6428,
                                                                       3398, 6458, 992, 1010,
                                                                       3578, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6788, 0, 3, 6458,
                                                                       3416, 6488, 1010, 1028,
                                                                       3614, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6848, 0, 3, 6488,
                                                                       3434, 6518, 1028, 1046,
                                                                       3650, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6908, 0, 3, 6518,
                                                                       3452, 6548, 1046, 1064,
                                                                       3686, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6968, 0, 3, 6548,
                                                                       3470, 6578, 1064, 1082,
                                                                       3722, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7028, 0, 3, 6608,
                                                                       3506, 6668, 1118, 1148,
                                                                       3758, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7128, 0, 3, 6668,
                                                                       3542, 6728, 1148, 1178,
                                                                       3818, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7228, 0, 3, 6728,
                                                                       3578, 6788, 1178, 1208,
                                                                       3878, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7328, 0, 3, 6788,
                                                                       3614, 6848, 1208, 1238,
                                                                       3938, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7428, 0, 3, 6848,
                                                                       3650, 6908, 1238, 1268,
                                                                       3998, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7528, 0, 3, 6908,
                                                                       3686, 6968, 1268, 1298,
                                                                       4058, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7628, 0, 3, 7028,
                                                                       3758, 7128, 1358, 1403,
                                                                       4118, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7778, 0, 3, 7128,
                                                                       3818, 7228, 1403, 1448,
                                                                       4208, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 7928, 0, 3, 7228,
                                                                       3878, 7328, 1448, 1493,
                                                                       4298, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8078, 0, 3, 7328,
                                                                       3938, 7428, 1493, 1538,
                                                                       4388, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8228, 0, 3, 7428,
                                                                       3998, 7528, 1538, 1583,
                                                                       4478, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 8378, 0, 3, 7628,
                                                                       4118, 7778, 1673, 1736,
                                                                       4568, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 8588, 0, 3, 7778,
                                                                       4208, 7928, 1736, 1799,
                                                                       4694, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 8798, 0, 3, 7928,
                                                                       4298, 8078, 1799, 1862,
                                                                       4820, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9008, 0, 3, 8078,
                                                                       4388, 8228, 1862, 1925,
                                                                       4946, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 9218, 0, 3, 8378,
                                                                       4568, 8588, 2051, 2135,
                                                                       5072, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 9498, 0, 3, 8588,
                                                                       4694, 8798, 2135, 2219,
                                                                       5240, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 9778, 0, 3, 8798,
                                                                       4820, 9008, 2219, 2303,
                                                                       5408, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 10058, 0, 3, 9218,
                                                                       5072, 9498, 2471, 2579,
                                                                       5576, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 10418, 0, 3, 9498,
                                                                       5240, 9778, 2579, 2687,
                                                                       5792, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 10778, 0, 3,
                                                                       10058, 5576, 10418, 2903,
                                                                       3038, 6008, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 11228, 9218, 280, ncols);

                    simdfunc::contract_primitives(buffer, 11704, 10058, 360, ncols);

                    simdfunc::contract_primitives(buffer, 12316, 10778, 450, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 11508, 11228, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 12064, 11704, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 12766, 12316, 45, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 13081, 11508, 12064, 7, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 13669, 12064, 12766, 7, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 14425, 13081, 13669, 7, nmax);

        simdtrf::transform_d_inner(buffer, 15601, 14425, 28, 7, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 15601, 35, nmax);
    }

    for (size_t m = 0; m < 455; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
