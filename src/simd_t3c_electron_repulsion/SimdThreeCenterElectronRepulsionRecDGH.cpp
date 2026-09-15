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


#include "SimdThreeCenterElectronRepulsionRecDGH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDG.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_dgh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_dgh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 20602, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 495 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 20602, 15782, 1740, dimensions);

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

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 103,
                                                                       113, 183, 198, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 309, 0, 3, 113,
                                                                       123, 198, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 330, 0, 3, 123,
                                                                       133, 213, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 351, 0, 3, 133,
                                                                       143, 228, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 372, 0, 3, 143,
                                                                       153, 243, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 393, 0, 3, 153,
                                                                       163, 258, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 414, 0, 3, 183,
                                                                       198, 288, 309, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 442, 0, 3, 198,
                                                                       213, 309, 330, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 470, 0, 3, 213,
                                                                       228, 330, 351, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 498, 0, 3, 228,
                                                                       243, 351, 372, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 526, 0, 3, 243,
                                                                       258, 372, 393, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 554, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 557, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 560, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 563, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 566, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 569, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 572, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 575, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 578, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 581, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 584, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 587, 3, 10, 25,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 596, 3, 11, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 605, 3, 12, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 614, 3, 13, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 623, 3, 14, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 632, 3, 15, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 641, 3, 16, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 650, 3, 17, 46,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 659, 3, 19, 49,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 677, 3, 22, 55,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 695, 3, 25, 61,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 713, 3, 28, 67,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 731, 3, 31, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 749, 3, 34, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 767, 3, 37, 85,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 785, 3, 40, 91,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 803, 3, 43, 97,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 821, 3, 49, 103,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 851, 3, 55, 113,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 881, 3, 61, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 911, 3, 67, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 941, 3, 73, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 971, 3, 79, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1001, 3, 85, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1031, 3, 91, 173,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1061, 3, 103, 183,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1106, 3, 113, 198,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1151, 3, 123, 213,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1196, 3, 133, 228,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1241, 3, 143, 243,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1286, 3, 153, 258,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1331, 3, 163, 273,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1376, 3, 183, 288,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1439, 3, 198, 309,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1502, 3, 213, 330,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1565, 3, 228, 351,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1628, 3, 243, 372,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1691, 3, 258, 393,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1754, 3, 288, 414,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1838, 3, 309, 442,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1922, 3, 330, 470,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2006, 3, 351, 498,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2090, 3, 372, 526,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2174, 3, 8, 9,
                                                                       560, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2180, 3, 9, 10,
                                                                       563, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2186, 3, 10, 11,
                                                                       566, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2192, 3, 11, 12,
                                                                       569, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2198, 3, 12, 13,
                                                                       572, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2204, 3, 13, 14,
                                                                       575, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2210, 3, 14, 15,
                                                                       578, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2216, 3, 15, 16,
                                                                       581, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2222, 3, 16, 17,
                                                                       584, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2228, 0, 3, 2174,
                                                                       560, 2180, 587, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2246, 0, 3, 2180,
                                                                       563, 2186, 596, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2264, 0, 3, 2186,
                                                                       566, 2192, 605, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2282, 0, 3, 2192,
                                                                       569, 2198, 614, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2300, 0, 3, 2198,
                                                                       572, 2204, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2318, 0, 3, 2204,
                                                                       575, 2210, 632, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2336, 0, 3, 2210,
                                                                       578, 2216, 641, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2354, 0, 3, 2216,
                                                                       581, 2222, 650, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 2228,
                                                                       587, 2246, 49, 55, 695,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2408, 0, 3, 2246,
                                                                       596, 2264, 55, 61, 713,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2444, 0, 3, 2264,
                                                                       605, 2282, 61, 67, 731,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2480, 0, 3, 2282,
                                                                       614, 2300, 67, 73, 749,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2516, 0, 3, 2300,
                                                                       623, 2318, 73, 79, 767,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2552, 0, 3, 2318,
                                                                       632, 2336, 79, 85, 785,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2588, 0, 3, 2336,
                                                                       641, 2354, 85, 91, 803,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2624, 0, 3, 2372,
                                                                       695, 2408, 103, 113, 881,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2684, 0, 3, 2408,
                                                                       713, 2444, 113, 123, 911,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2744, 0, 3, 2444,
                                                                       731, 2480, 123, 133, 941,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2804, 0, 3, 2480,
                                                                       749, 2516, 133, 143, 971,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2864, 0, 3, 2516,
                                                                       767, 2552, 143, 153, 1001,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2924, 0, 3, 2552,
                                                                       785, 2588, 153, 163, 1031,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2984, 0, 3, 2624,
                                                                       881, 2684, 183, 198, 1151,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3074, 0, 3, 2684,
                                                                       911, 2744, 198, 213, 1196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3164, 0, 3, 2744,
                                                                       941, 2804, 213, 228, 1241,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3254, 0, 3, 2804,
                                                                       971, 2864, 228, 243, 1286,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3344, 0, 3, 2864,
                                                                       1001, 2924, 243, 258,
                                                                       1331, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3434, 0, 3, 2984,
                                                                       1151, 3074, 288, 309,
                                                                       1502, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3560, 0, 3, 3074,
                                                                       1196, 3164, 309, 330,
                                                                       1565, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3686, 0, 3, 3164,
                                                                       1241, 3254, 330, 351,
                                                                       1628, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3812, 0, 3, 3254,
                                                                       1286, 3344, 351, 372,
                                                                       1691, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 3938, 0, 3, 3434,
                                                                       1502, 3560, 414, 442,
                                                                       1922, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 4106, 0, 3, 3560,
                                                                       1565, 3686, 442, 470,
                                                                       2006, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 4274, 0, 3, 3686,
                                                                       1628, 3812, 470, 498,
                                                                       2090, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4442, 3, 554, 557,
                                                                       2174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4452, 3, 557, 560,
                                                                       2180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4462, 3, 560, 563,
                                                                       2186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4472, 3, 563, 566,
                                                                       2192, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4482, 3, 566, 569,
                                                                       2198, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4492, 3, 569, 572,
                                                                       2204, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4502, 3, 572, 575,
                                                                       2210, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4512, 3, 575, 578,
                                                                       2216, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4522, 3, 578, 581,
                                                                       2222, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4532, 0, 3, 4442,
                                                                       2174, 4452, 2228, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4562, 0, 3, 4452,
                                                                       2180, 4462, 2246, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4592, 0, 3, 4462,
                                                                       2186, 4472, 2264, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4622, 0, 3, 4472,
                                                                       2192, 4482, 2282, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4652, 0, 3, 4482,
                                                                       2198, 4492, 2300, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4682, 0, 3, 4492,
                                                                       2204, 4502, 2318, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4712, 0, 3, 4502,
                                                                       2210, 4512, 2336, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4742, 0, 3, 4512,
                                                                       2216, 4522, 2354, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4772, 0, 3, 4532,
                                                                       2228, 4562, 659, 677,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4832, 0, 3, 4562,
                                                                       2246, 4592, 677, 695,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4892, 0, 3, 4592,
                                                                       2264, 4622, 695, 713,
                                                                       2444, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4952, 0, 3, 4622,
                                                                       2282, 4652, 713, 731,
                                                                       2480, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5012, 0, 3, 4652,
                                                                       2300, 4682, 731, 749,
                                                                       2516, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5072, 0, 3, 4682,
                                                                       2318, 4712, 749, 767,
                                                                       2552, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5132, 0, 3, 4712,
                                                                       2336, 4742, 767, 785,
                                                                       2588, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5192, 0, 3, 4772,
                                                                       2372, 4832, 821, 851,
                                                                       2624, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5292, 0, 3, 4832,
                                                                       2408, 4892, 851, 881,
                                                                       2684, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5392, 0, 3, 4892,
                                                                       2444, 4952, 881, 911,
                                                                       2744, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5492, 0, 3, 4952,
                                                                       2480, 5012, 911, 941,
                                                                       2804, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5592, 0, 3, 5012,
                                                                       2516, 5072, 941, 971,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5692, 0, 3, 5072,
                                                                       2552, 5132, 971, 1001,
                                                                       2924, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5792, 0, 3, 5192,
                                                                       2624, 5292, 1061, 1106,
                                                                       2984, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5942, 0, 3, 5292,
                                                                       2684, 5392, 1106, 1151,
                                                                       3074, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6092, 0, 3, 5392,
                                                                       2744, 5492, 1151, 1196,
                                                                       3164, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6242, 0, 3, 5492,
                                                                       2804, 5592, 1196, 1241,
                                                                       3254, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6392, 0, 3, 5592,
                                                                       2864, 5692, 1241, 1286,
                                                                       3344, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 6542, 0, 3, 5792,
                                                                       2984, 5942, 1376, 1439,
                                                                       3434, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 6752, 0, 3, 5942,
                                                                       3074, 6092, 1439, 1502,
                                                                       3560, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 6962, 0, 3, 6092,
                                                                       3164, 6242, 1502, 1565,
                                                                       3686, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7172, 0, 3, 6242,
                                                                       3254, 6392, 1565, 1628,
                                                                       3812, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 7382, 0, 3, 6542,
                                                                       3434, 6752, 1754, 1838,
                                                                       3938, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 7662, 0, 3, 6752,
                                                                       3560, 6962, 1838, 1922,
                                                                       4106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 7942, 0, 3, 6962,
                                                                       3686, 7172, 1922, 2006,
                                                                       4274, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8222, 3, 2174,
                                                                       2180, 4462, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8237, 3, 2180,
                                                                       2186, 4472, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8252, 3, 2186,
                                                                       2192, 4482, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8267, 3, 2192,
                                                                       2198, 4492, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8282, 3, 2198,
                                                                       2204, 4502, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8297, 3, 2204,
                                                                       2210, 4512, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8312, 3, 2210,
                                                                       2216, 4522, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8327, 0, 3, 8222,
                                                                       4462, 8237, 2228, 2246,
                                                                       4592, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8372, 0, 3, 8237,
                                                                       4472, 8252, 2246, 2264,
                                                                       4622, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8417, 0, 3, 8252,
                                                                       4482, 8267, 2264, 2282,
                                                                       4652, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8462, 0, 3, 8267,
                                                                       4492, 8282, 2282, 2300,
                                                                       4682, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8507, 0, 3, 8282,
                                                                       4502, 8297, 2300, 2318,
                                                                       4712, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8552, 0, 3, 8297,
                                                                       4512, 8312, 2318, 2336,
                                                                       4742, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8597, 0, 3, 8327,
                                                                       4592, 8372, 2372, 2408,
                                                                       4892, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8687, 0, 3, 8372,
                                                                       4622, 8417, 2408, 2444,
                                                                       4952, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8777, 0, 3, 8417,
                                                                       4652, 8462, 2444, 2480,
                                                                       5012, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8867, 0, 3, 8462,
                                                                       4682, 8507, 2480, 2516,
                                                                       5072, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8957, 0, 3, 8507,
                                                                       4712, 8552, 2516, 2552,
                                                                       5132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9047, 0, 3, 8597,
                                                                       4892, 8687, 2624, 2684,
                                                                       5392, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9197, 0, 3, 8687,
                                                                       4952, 8777, 2684, 2744,
                                                                       5492, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9347, 0, 3, 8777,
                                                                       5012, 8867, 2744, 2804,
                                                                       5592, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9497, 0, 3, 8867,
                                                                       5072, 8957, 2804, 2864,
                                                                       5692, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 9647, 0, 3, 9047,
                                                                       5392, 9197, 2984, 3074,
                                                                       6092, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 9872, 0, 3, 9197,
                                                                       5492, 9347, 3074, 3164,
                                                                       6242, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 10097, 0, 3, 9347,
                                                                       5592, 9497, 3164, 3254,
                                                                       6392, ncols, gamma, p,
                                                                       q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 10322, 0, 3, 9647,
                                                                       6092, 9872, 3434, 3560,
                                                                       6962, ncols, gamma, p,
                                                                       q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 10637, 0, 3, 9872,
                                                                       6242, 10097, 3560, 3686,
                                                                       7172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 10952, 0, 3,
                                                                       10322, 6962, 10637, 3938,
                                                                       4106, 7942, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11372, 3, 4442,
                                                                       4452, 8222, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11393, 3, 4452,
                                                                       4462, 8237, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11414, 3, 4462,
                                                                       4472, 8252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11435, 3, 4472,
                                                                       4482, 8267, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11456, 3, 4482,
                                                                       4492, 8282, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11477, 3, 4492,
                                                                       4502, 8297, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11498, 3, 4502,
                                                                       4512, 8312, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11519, 0, 3,
                                                                       11372, 8222, 11393, 4532,
                                                                       4562, 8327, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11582, 0, 3,
                                                                       11393, 8237, 11414, 4562,
                                                                       4592, 8372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11645, 0, 3,
                                                                       11414, 8252, 11435, 4592,
                                                                       4622, 8417, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11708, 0, 3,
                                                                       11435, 8267, 11456, 4622,
                                                                       4652, 8462, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11771, 0, 3,
                                                                       11456, 8282, 11477, 4652,
                                                                       4682, 8507, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11834, 0, 3,
                                                                       11477, 8297, 11498, 4682,
                                                                       4712, 8552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 11897, 0, 3,
                                                                       11519, 8327, 11582, 4772,
                                                                       4832, 8597, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 12023, 0, 3,
                                                                       11582, 8372, 11645, 4832,
                                                                       4892, 8687, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 12149, 0, 3,
                                                                       11645, 8417, 11708, 4892,
                                                                       4952, 8777, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 12275, 0, 3,
                                                                       11708, 8462, 11771, 4952,
                                                                       5012, 8867, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 12401, 0, 3,
                                                                       11771, 8507, 11834, 5012,
                                                                       5072, 8957, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 12527, 0, 3,
                                                                       11897, 8597, 12023, 5192,
                                                                       5292, 9047, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 12737, 0, 3,
                                                                       12023, 8687, 12149, 5292,
                                                                       5392, 9197, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 12947, 0, 3,
                                                                       12149, 8777, 12275, 5392,
                                                                       5492, 9347, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 13157, 0, 3,
                                                                       12275, 8867, 12401, 5492,
                                                                       5592, 9497, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 13367, 0, 3,
                                                                       12527, 9047, 12737, 5792,
                                                                       5942, 9647, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 13682, 0, 3,
                                                                       12737, 9197, 12947, 5942,
                                                                       6092, 9872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 13997, 0, 3,
                                                                       12947, 9347, 13157, 6092,
                                                                       6242, 10097, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 14312, 0, 3,
                                                                       13367, 9647, 13682, 6542,
                                                                       6752, 10322, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 14753, 0, 3,
                                                                       13682, 9872, 13997, 6752,
                                                                       6962, 10637, ncols, gamma,
                                                                       p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 15194, 0, 3,
                                                                       14312, 10322, 14753, 7382,
                                                                       7662, 10952, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 15782, 13367, 315, ncols);

                    simdfunc::contract_primitives(buffer, 16262, 14312, 441, ncols);

                    simdfunc::contract_primitives(buffer, 16934, 15194, 588, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 16097, 15782, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 16703, 16262, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 17522, 16934, 28, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 17830, 16097, 16703, 11, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 18325, 16703, 17522, 11, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 19018, 17830, 18325, 11, nmax);

        simdtrf::transform_g_inner(buffer, 20008, 19018, 6, 11, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 20008, 99, nmax);
    }

    for (size_t m = 0; m < 495; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
