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


#include "SimdThreeCenterElectronRepulsionRecSGK.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_sgk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_sgk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 14755, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 135 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 14755, 13990, 540, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

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

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 19, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 49, 0, 3, 8, 9,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 55, 0, 3, 9, 10,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 61, 0, 3, 10, 11,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 67, 0, 3, 11, 12,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 73, 0, 3, 12, 13,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 79, 0, 3, 13, 14,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 85, 0, 3, 14, 15,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 91, 0, 3, 15, 16,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 97, 0, 3, 16, 17,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 103, 0, 3, 19, 22,
                                                                       49, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 113, 0, 3, 22, 25,
                                                                       55, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 123, 0, 3, 25, 28,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 133, 0, 3, 28, 31,
                                                                       67, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 143, 0, 3, 31, 34,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 34, 37,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 37, 40,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 40, 43,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 49, 55,
                                                                       103, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 198, 0, 3, 55, 61,
                                                                       113, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 61, 67,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 67, 73,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 73, 79,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 79, 85,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 85, 91,
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

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 321, 3, 10, 25,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 330, 3, 11, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 339, 3, 12, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 348, 3, 13, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 357, 3, 14, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 366, 3, 15, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 375, 3, 16, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 384, 3, 17, 46,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 393, 3, 19, 49,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 411, 3, 22, 55,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 429, 3, 25, 61,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 447, 3, 28, 67,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 465, 3, 31, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 483, 3, 34, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 501, 3, 37, 85,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 519, 3, 40, 91,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 537, 3, 43, 97,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 555, 3, 49, 103,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 585, 3, 55, 113,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 615, 3, 61, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 645, 3, 67, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 675, 3, 73, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 705, 3, 79, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 735, 3, 85, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 765, 3, 91, 173,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 795, 3, 103, 183,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 840, 3, 113, 198,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 885, 3, 123, 213,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 930, 3, 133, 228,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 975, 3, 143, 243,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1020, 3, 153, 258,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1065, 3, 163, 273,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1110, 3, 8, 9,
                                                                       294, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1116, 3, 9, 10,
                                                                       297, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1122, 3, 10, 11,
                                                                       300, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1128, 3, 11, 12,
                                                                       303, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1134, 3, 12, 13,
                                                                       306, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1140, 3, 13, 14,
                                                                       309, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1146, 3, 14, 15,
                                                                       312, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1152, 3, 15, 16,
                                                                       315, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1158, 3, 16, 17,
                                                                       318, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1164, 0, 3, 1110,
                                                                       294, 1116, 321, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1182, 0, 3, 1116,
                                                                       297, 1122, 330, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1200, 0, 3, 1122,
                                                                       300, 1128, 339, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1218, 0, 3, 1128,
                                                                       303, 1134, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1236, 0, 3, 1134,
                                                                       306, 1140, 357, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1254, 0, 3, 1140,
                                                                       309, 1146, 366, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1272, 0, 3, 1146,
                                                                       312, 1152, 375, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1290, 0, 3, 1152,
                                                                       315, 1158, 384, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1308, 0, 3, 1164,
                                                                       321, 1182, 49, 55, 429,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1344, 0, 3, 1182,
                                                                       330, 1200, 55, 61, 447,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 1200,
                                                                       339, 1218, 61, 67, 465,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1416, 0, 3, 1218,
                                                                       348, 1236, 67, 73, 483,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1452, 0, 3, 1236,
                                                                       357, 1254, 73, 79, 501,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1488, 0, 3, 1254,
                                                                       366, 1272, 79, 85, 519,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1524, 0, 3, 1272,
                                                                       375, 1290, 85, 91, 537,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1560, 0, 3, 1308,
                                                                       429, 1344, 103, 113, 615,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1620, 0, 3, 1344,
                                                                       447, 1380, 113, 123, 645,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1680, 0, 3, 1380,
                                                                       465, 1416, 123, 133, 675,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1740, 0, 3, 1416,
                                                                       483, 1452, 133, 143, 705,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 1452,
                                                                       501, 1488, 143, 153, 735,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1860, 0, 3, 1488,
                                                                       519, 1524, 153, 163, 765,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 1920, 0, 3, 1560,
                                                                       615, 1620, 183, 198, 885,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2010, 0, 3, 1620,
                                                                       645, 1680, 198, 213, 930,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2100, 0, 3, 1680,
                                                                       675, 1740, 213, 228, 975,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2190, 0, 3, 1740,
                                                                       705, 1800, 228, 243, 1020,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2280, 0, 3, 1800,
                                                                       735, 1860, 243, 258, 1065,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2370, 3, 288, 291,
                                                                       1110, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2380, 3, 291, 294,
                                                                       1116, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2390, 3, 294, 297,
                                                                       1122, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2400, 3, 297, 300,
                                                                       1128, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2410, 3, 300, 303,
                                                                       1134, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2420, 3, 303, 306,
                                                                       1140, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2430, 3, 306, 309,
                                                                       1146, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2440, 3, 309, 312,
                                                                       1152, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2450, 3, 312, 315,
                                                                       1158, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2460, 0, 3, 2370,
                                                                       1110, 2380, 1164, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2490, 0, 3, 2380,
                                                                       1116, 2390, 1182, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2520, 0, 3, 2390,
                                                                       1122, 2400, 1200, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2550, 0, 3, 2400,
                                                                       1128, 2410, 1218, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2580, 0, 3, 2410,
                                                                       1134, 2420, 1236, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2610, 0, 3, 2420,
                                                                       1140, 2430, 1254, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2640, 0, 3, 2430,
                                                                       1146, 2440, 1272, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2670, 0, 3, 2440,
                                                                       1152, 2450, 1290, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2700, 0, 3, 2460,
                                                                       1164, 2490, 393, 411,
                                                                       1308, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2760, 0, 3, 2490,
                                                                       1182, 2520, 411, 429,
                                                                       1344, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2820, 0, 3, 2520,
                                                                       1200, 2550, 429, 447,
                                                                       1380, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2880, 0, 3, 2550,
                                                                       1218, 2580, 447, 465,
                                                                       1416, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2940, 0, 3, 2580,
                                                                       1236, 2610, 465, 483,
                                                                       1452, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3000, 0, 3, 2610,
                                                                       1254, 2640, 483, 501,
                                                                       1488, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3060, 0, 3, 2640,
                                                                       1272, 2670, 501, 519,
                                                                       1524, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3120, 0, 3, 2700,
                                                                       1308, 2760, 555, 585,
                                                                       1560, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3220, 0, 3, 2760,
                                                                       1344, 2820, 585, 615,
                                                                       1620, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3320, 0, 3, 2820,
                                                                       1380, 2880, 615, 645,
                                                                       1680, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3420, 0, 3, 2880,
                                                                       1416, 2940, 645, 675,
                                                                       1740, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3520, 0, 3, 2940,
                                                                       1452, 3000, 675, 705,
                                                                       1800, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3620, 0, 3, 3000,
                                                                       1488, 3060, 705, 735,
                                                                       1860, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 3720, 0, 3, 3120,
                                                                       1560, 3220, 795, 840,
                                                                       1920, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 3870, 0, 3, 3220,
                                                                       1620, 3320, 840, 885,
                                                                       2010, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 4020, 0, 3, 3320,
                                                                       1680, 3420, 885, 930,
                                                                       2100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 4170, 0, 3, 3420,
                                                                       1740, 3520, 930, 975,
                                                                       2190, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 4320, 0, 3, 3520,
                                                                       1800, 3620, 975, 1020,
                                                                       2280, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4470, 3, 1110,
                                                                       1116, 2390, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4485, 3, 1116,
                                                                       1122, 2400, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4500, 3, 1122,
                                                                       1128, 2410, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4515, 3, 1128,
                                                                       1134, 2420, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4530, 3, 1134,
                                                                       1140, 2430, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4545, 3, 1140,
                                                                       1146, 2440, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4560, 3, 1146,
                                                                       1152, 2450, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4575, 0, 3, 4470,
                                                                       2390, 4485, 1164, 1182,
                                                                       2520, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4620, 0, 3, 4485,
                                                                       2400, 4500, 1182, 1200,
                                                                       2550, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4665, 0, 3, 4500,
                                                                       2410, 4515, 1200, 1218,
                                                                       2580, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4710, 0, 3, 4515,
                                                                       2420, 4530, 1218, 1236,
                                                                       2610, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4755, 0, 3, 4530,
                                                                       2430, 4545, 1236, 1254,
                                                                       2640, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4800, 0, 3, 4545,
                                                                       2440, 4560, 1254, 1272,
                                                                       2670, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 4845, 0, 3, 4575,
                                                                       2520, 4620, 1308, 1344,
                                                                       2820, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 4935, 0, 3, 4620,
                                                                       2550, 4665, 1344, 1380,
                                                                       2880, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5025, 0, 3, 4665,
                                                                       2580, 4710, 1380, 1416,
                                                                       2940, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5115, 0, 3, 4710,
                                                                       2610, 4755, 1416, 1452,
                                                                       3000, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5205, 0, 3, 4755,
                                                                       2640, 4800, 1452, 1488,
                                                                       3060, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 5295, 0, 3, 4845,
                                                                       2820, 4935, 1560, 1620,
                                                                       3320, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 5445, 0, 3, 4935,
                                                                       2880, 5025, 1620, 1680,
                                                                       3420, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 5595, 0, 3, 5025,
                                                                       2940, 5115, 1680, 1740,
                                                                       3520, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 5745, 0, 3, 5115,
                                                                       3000, 5205, 1740, 1800,
                                                                       3620, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 5895, 0, 3, 5295,
                                                                       3320, 5445, 1920, 2010,
                                                                       4020, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 6120, 0, 3, 5445,
                                                                       3420, 5595, 2010, 2100,
                                                                       4170, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 6345, 0, 3, 5595,
                                                                       3520, 5745, 2100, 2190,
                                                                       4320, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6570, 3, 2370,
                                                                       2380, 4470, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6591, 3, 2380,
                                                                       2390, 4485, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6612, 3, 2390,
                                                                       2400, 4500, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6633, 3, 2400,
                                                                       2410, 4515, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6654, 3, 2410,
                                                                       2420, 4530, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6675, 3, 2420,
                                                                       2430, 4545, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6696, 3, 2430,
                                                                       2440, 4560, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 6717, 0, 3, 6570,
                                                                       4470, 6591, 2460, 2490,
                                                                       4575, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 6780, 0, 3, 6591,
                                                                       4485, 6612, 2490, 2520,
                                                                       4620, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 6843, 0, 3, 6612,
                                                                       4500, 6633, 2520, 2550,
                                                                       4665, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 6906, 0, 3, 6633,
                                                                       4515, 6654, 2550, 2580,
                                                                       4710, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 6969, 0, 3, 6654,
                                                                       4530, 6675, 2580, 2610,
                                                                       4755, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7032, 0, 3, 6675,
                                                                       4545, 6696, 2610, 2640,
                                                                       4800, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 7095, 0, 3, 6717,
                                                                       4575, 6780, 2700, 2760,
                                                                       4845, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 7221, 0, 3, 6780,
                                                                       4620, 6843, 2760, 2820,
                                                                       4935, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 7347, 0, 3, 6843,
                                                                       4665, 6906, 2820, 2880,
                                                                       5025, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 7473, 0, 3, 6906,
                                                                       4710, 6969, 2880, 2940,
                                                                       5115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 7599, 0, 3, 6969,
                                                                       4755, 7032, 2940, 3000,
                                                                       5205, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 7725, 0, 3, 7095,
                                                                       4845, 7221, 3120, 3220,
                                                                       5295, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 7935, 0, 3, 7221,
                                                                       4935, 7347, 3220, 3320,
                                                                       5445, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 8145, 0, 3, 7347,
                                                                       5025, 7473, 3320, 3420,
                                                                       5595, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 8355, 0, 3, 7473,
                                                                       5115, 7599, 3420, 3520,
                                                                       5745, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 8565, 0, 3, 7725,
                                                                       5295, 7935, 3720, 3870,
                                                                       5895, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 8880, 0, 3, 7935,
                                                                       5445, 8145, 3870, 4020,
                                                                       6120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 9195, 0, 3, 8145,
                                                                       5595, 8355, 4020, 4170,
                                                                       6345, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 9510, 3, 4470,
                                                                       4485, 6612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 9538, 3, 4485,
                                                                       4500, 6633, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 9566, 3, 4500,
                                                                       4515, 6654, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 9594, 3, 4515,
                                                                       4530, 6675, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 9622, 3, 4530,
                                                                       4545, 6696, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 9650, 0, 3, 9510,
                                                                       6612, 9538, 4575, 4620,
                                                                       6843, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 9734, 0, 3, 9538,
                                                                       6633, 9566, 4620, 4665,
                                                                       6906, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 9818, 0, 3, 9566,
                                                                       6654, 9594, 4665, 4710,
                                                                       6969, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 9902, 0, 3, 9594,
                                                                       6675, 9622, 4710, 4755,
                                                                       7032, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 9986, 0, 3, 9650,
                                                                       6843, 9734, 4845, 4935,
                                                                       7347, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 10154, 0, 3, 9734,
                                                                       6906, 9818, 4935, 5025,
                                                                       7473, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 10322, 0, 3, 9818,
                                                                       6969, 9902, 5025, 5115,
                                                                       7599, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 10490, 0, 3, 9986,
                                                                       7347, 10154, 5295, 5445,
                                                                       8145, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 10770, 0, 3,
                                                                       10154, 7473, 10322, 5445,
                                                                       5595, 8355, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 11050, 0, 3,
                                                                       10490, 8145, 10770, 5895,
                                                                       6120, 9195, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 11470, 3, 6570,
                                                                       6591, 9510, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 11506, 3, 6591,
                                                                       6612, 9538, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 11542, 3, 6612,
                                                                       6633, 9566, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 11578, 3, 6633,
                                                                       6654, 9594, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 11614, 3, 6654,
                                                                       6675, 9622, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 11650, 0, 3,
                                                                       11470, 9510, 11506, 6717,
                                                                       6780, 9650, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 11758, 0, 3,
                                                                       11506, 9538, 11542, 6780,
                                                                       6843, 9734, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 11866, 0, 3,
                                                                       11542, 9566, 11578, 6843,
                                                                       6906, 9818, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 11974, 0, 3,
                                                                       11578, 9594, 11614, 6906,
                                                                       6969, 9902, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 12082, 0, 3,
                                                                       11650, 9650, 11758, 7095,
                                                                       7221, 9986, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 12298, 0, 3,
                                                                       11758, 9734, 11866, 7221,
                                                                       7347, 10154, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 12514, 0, 3,
                                                                       11866, 9818, 11974, 7347,
                                                                       7473, 10322, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 12730, 0, 3,
                                                                       12082, 9986, 12298, 7725,
                                                                       7935, 10490, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 13090, 0, 3,
                                                                       12298, 10154, 12514, 7935,
                                                                       8145, 10770, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 13450, 0, 3,
                                                                       12730, 10490, 13090, 8565,
                                                                       8880, 11050, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 13990, 13450, 540, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 14530, 13990, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 14530, 15, nmax);
    }

    for (size_t m = 0; m < 135; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
