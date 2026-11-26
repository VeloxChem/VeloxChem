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


#include "SimdThreeCenterElectronRepulsionRsRecPDK.hpp"

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
#include "SimdTransferPD.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_pdk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_pdk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 17169, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 17169, 14772, 1482, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

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

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 73, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 76, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 79, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 82, 0, 3, 7, 8,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 88, 0, 3, 8, 9,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 94, 0, 3, 9, 10,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 100, 0, 3, 10, 11,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 106, 0, 3, 11, 12,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 112, 0, 3, 12, 13,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 118, 0, 3, 13, 14,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 124, 0, 3, 14, 15,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 130, 0, 3, 18, 19,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 136, 0, 3, 19, 20,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 142, 0, 3, 20, 21,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 148, 0, 3, 21, 22,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 154, 0, 3, 22, 23,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 160, 0, 3, 23, 24,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 166, 0, 3, 24, 25,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 172, 0, 3, 25, 26,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 178, 0, 3, 28, 31,
                                                                       82, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 188, 0, 3, 31, 34,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 198, 0, 3, 34, 37,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 208, 0, 3, 37, 40,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 40, 43,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 43, 46,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 46, 49,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 248, 0, 3, 55, 58,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 58, 61,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 61, 64,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 64, 67,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 67, 70,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 70, 73,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 73, 76,
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

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 378, 3, 9, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 387, 3, 10, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 396, 3, 11, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 405, 3, 12, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 414, 3, 13, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 423, 3, 14, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 432, 3, 15, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 441, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 450, 3, 21, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 459, 3, 22, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 468, 3, 23, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 477, 3, 24, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 486, 3, 25, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 495, 3, 26, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 504, 3, 28, 82,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 522, 3, 31, 88,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 540, 3, 34, 94,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 558, 3, 37, 100,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 576, 3, 40, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 594, 3, 43, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 612, 3, 46, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 630, 3, 49, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 648, 3, 55, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 666, 3, 58, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 684, 3, 61, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 702, 3, 64, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 720, 3, 67, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 738, 3, 70, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 756, 3, 73, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 774, 3, 76, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 792, 3, 82, 178,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 822, 3, 88, 188,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 852, 3, 94, 198,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 882, 3, 100, 208,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 912, 3, 106, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 942, 3, 112, 228,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 972, 3, 118, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1002, 3, 130, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1032, 3, 136, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1062, 3, 142, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1092, 3, 148, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1122, 3, 154, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1152, 3, 160, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1182, 3, 166, 308,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1212, 3, 7, 8,
                                                                       324, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1218, 3, 8, 9,
                                                                       327, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1224, 3, 9, 10,
                                                                       330, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1230, 3, 10, 11,
                                                                       333, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1236, 3, 11, 12,
                                                                       336, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1242, 3, 12, 13,
                                                                       339, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1248, 3, 13, 14,
                                                                       342, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1254, 3, 14, 15,
                                                                       345, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1260, 3, 18, 19,
                                                                       354, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1266, 3, 19, 20,
                                                                       357, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1272, 3, 20, 21,
                                                                       360, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1278, 3, 21, 22,
                                                                       363, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1284, 3, 22, 23,
                                                                       366, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1290, 3, 23, 24,
                                                                       369, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1296, 3, 24, 25,
                                                                       372, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1302, 3, 25, 26,
                                                                       375, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1308, 0, 3, 1212,
                                                                       324, 1218, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1326, 0, 3, 1218,
                                                                       327, 1224, 387, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1344, 0, 3, 1224,
                                                                       330, 1230, 396, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1362, 0, 3, 1230,
                                                                       333, 1236, 405, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 1236,
                                                                       336, 1242, 414, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1398, 0, 3, 1242,
                                                                       339, 1248, 423, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1416, 0, 3, 1248,
                                                                       342, 1254, 432, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1434, 0, 3, 1260,
                                                                       354, 1266, 441, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1452, 0, 3, 1266,
                                                                       357, 1272, 450, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1470, 0, 3, 1272,
                                                                       360, 1278, 459, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1488, 0, 3, 1278,
                                                                       363, 1284, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1506, 0, 3, 1284,
                                                                       366, 1290, 477, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1524, 0, 3, 1290,
                                                                       369, 1296, 486, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1542, 0, 3, 1296,
                                                                       372, 1302, 495, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1560, 0, 3, 1308,
                                                                       378, 1326, 82, 88, 540,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1596, 0, 3, 1326,
                                                                       387, 1344, 88, 94, 558,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 1344,
                                                                       396, 1362, 94, 100, 576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1668, 0, 3, 1362,
                                                                       405, 1380, 100, 106, 594,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1704, 0, 3, 1380,
                                                                       414, 1398, 106, 112, 612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1740, 0, 3, 1398,
                                                                       423, 1416, 112, 118, 630,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1776, 0, 3, 1434,
                                                                       441, 1452, 130, 136, 684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1812, 0, 3, 1452,
                                                                       450, 1470, 136, 142, 702,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1848, 0, 3, 1470,
                                                                       459, 1488, 142, 148, 720,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 1488,
                                                                       468, 1506, 148, 154, 738,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1920, 0, 3, 1506,
                                                                       477, 1524, 154, 160, 756,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1956, 0, 3, 1524,
                                                                       486, 1542, 160, 166, 774,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1992, 0, 3, 1560,
                                                                       540, 1596, 178, 188, 852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2052, 0, 3, 1596,
                                                                       558, 1632, 188, 198, 882,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2112, 0, 3, 1632,
                                                                       576, 1668, 198, 208, 912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2172, 0, 3, 1668,
                                                                       594, 1704, 208, 218, 942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2232, 0, 3, 1704,
                                                                       612, 1740, 218, 228, 972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2292, 0, 3, 1776,
                                                                       684, 1812, 248, 258, 1062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2352, 0, 3, 1812,
                                                                       702, 1848, 258, 268, 1092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2412, 0, 3, 1848,
                                                                       720, 1884, 268, 278, 1122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2472, 0, 3, 1884,
                                                                       738, 1920, 278, 288, 1152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2532, 0, 3, 1920,
                                                                       756, 1956, 288, 298, 1182,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2592, 3, 318, 321,
                                                                       1212, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2602, 3, 321, 324,
                                                                       1218, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2612, 3, 324, 327,
                                                                       1224, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2622, 3, 327, 330,
                                                                       1230, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2632, 3, 330, 333,
                                                                       1236, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2642, 3, 333, 336,
                                                                       1242, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2652, 3, 336, 339,
                                                                       1248, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2662, 3, 339, 342,
                                                                       1254, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2672, 3, 348, 351,
                                                                       1260, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2682, 3, 351, 354,
                                                                       1266, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2692, 3, 354, 357,
                                                                       1272, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2702, 3, 357, 360,
                                                                       1278, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2712, 3, 360, 363,
                                                                       1284, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2722, 3, 363, 366,
                                                                       1290, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2732, 3, 366, 369,
                                                                       1296, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2742, 3, 369, 372,
                                                                       1302, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2752, 0, 3, 2592,
                                                                       1212, 2602, 1308, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2782, 0, 3, 2602,
                                                                       1218, 2612, 1326, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2812, 0, 3, 2612,
                                                                       1224, 2622, 1344, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2842, 0, 3, 2622,
                                                                       1230, 2632, 1362, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2872, 0, 3, 2632,
                                                                       1236, 2642, 1380, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2902, 0, 3, 2642,
                                                                       1242, 2652, 1398, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2932, 0, 3, 2652,
                                                                       1248, 2662, 1416, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2962, 0, 3, 2672,
                                                                       1260, 2682, 1434, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2992, 0, 3, 2682,
                                                                       1266, 2692, 1452, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3022, 0, 3, 2692,
                                                                       1272, 2702, 1470, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3052, 0, 3, 2702,
                                                                       1278, 2712, 1488, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3082, 0, 3, 2712,
                                                                       1284, 2722, 1506, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3112, 0, 3, 2722,
                                                                       1290, 2732, 1524, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3142, 0, 3, 2732,
                                                                       1296, 2742, 1542, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3172, 0, 3, 2752,
                                                                       1308, 2782, 504, 522,
                                                                       1560, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3232, 0, 3, 2782,
                                                                       1326, 2812, 522, 540,
                                                                       1596, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3292, 0, 3, 2812,
                                                                       1344, 2842, 540, 558,
                                                                       1632, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3352, 0, 3, 2842,
                                                                       1362, 2872, 558, 576,
                                                                       1668, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3412, 0, 3, 2872,
                                                                       1380, 2902, 576, 594,
                                                                       1704, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3472, 0, 3, 2902,
                                                                       1398, 2932, 594, 612,
                                                                       1740, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3532, 0, 3, 2962,
                                                                       1434, 2992, 648, 666,
                                                                       1776, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3592, 0, 3, 2992,
                                                                       1452, 3022, 666, 684,
                                                                       1812, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3652, 0, 3, 3022,
                                                                       1470, 3052, 684, 702,
                                                                       1848, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3712, 0, 3, 3052,
                                                                       1488, 3082, 702, 720,
                                                                       1884, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3772, 0, 3, 3082,
                                                                       1506, 3112, 720, 738,
                                                                       1920, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3832, 0, 3, 3112,
                                                                       1524, 3142, 738, 756,
                                                                       1956, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3892, 0, 3, 3172,
                                                                       1560, 3232, 792, 822,
                                                                       1992, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3992, 0, 3, 3232,
                                                                       1596, 3292, 822, 852,
                                                                       2052, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4092, 0, 3, 3292,
                                                                       1632, 3352, 852, 882,
                                                                       2112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4192, 0, 3, 3352,
                                                                       1668, 3412, 882, 912,
                                                                       2172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4292, 0, 3, 3412,
                                                                       1704, 3472, 912, 942,
                                                                       2232, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4392, 0, 3, 3532,
                                                                       1776, 3592, 1002, 1032,
                                                                       2292, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4492, 0, 3, 3592,
                                                                       1812, 3652, 1032, 1062,
                                                                       2352, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4592, 0, 3, 3652,
                                                                       1848, 3712, 1062, 1092,
                                                                       2412, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4692, 0, 3, 3712,
                                                                       1884, 3772, 1092, 1122,
                                                                       2472, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4792, 0, 3, 3772,
                                                                       1920, 3832, 1122, 1152,
                                                                       2532, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4892, 3, 1212,
                                                                       1218, 2612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4907, 3, 1218,
                                                                       1224, 2622, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4922, 3, 1224,
                                                                       1230, 2632, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4937, 3, 1230,
                                                                       1236, 2642, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4952, 3, 1236,
                                                                       1242, 2652, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4967, 3, 1242,
                                                                       1248, 2662, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4982, 3, 1260,
                                                                       1266, 2692, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4997, 3, 1266,
                                                                       1272, 2702, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5012, 3, 1272,
                                                                       1278, 2712, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5027, 3, 1278,
                                                                       1284, 2722, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5042, 3, 1284,
                                                                       1290, 2732, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5057, 3, 1290,
                                                                       1296, 2742, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5072, 0, 3, 4892,
                                                                       2612, 4907, 1308, 1326,
                                                                       2812, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5117, 0, 3, 4907,
                                                                       2622, 4922, 1326, 1344,
                                                                       2842, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5162, 0, 3, 4922,
                                                                       2632, 4937, 1344, 1362,
                                                                       2872, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5207, 0, 3, 4937,
                                                                       2642, 4952, 1362, 1380,
                                                                       2902, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5252, 0, 3, 4952,
                                                                       2652, 4967, 1380, 1398,
                                                                       2932, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5297, 0, 3, 4982,
                                                                       2692, 4997, 1434, 1452,
                                                                       3022, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5342, 0, 3, 4997,
                                                                       2702, 5012, 1452, 1470,
                                                                       3052, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5387, 0, 3, 5012,
                                                                       2712, 5027, 1470, 1488,
                                                                       3082, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5432, 0, 3, 5027,
                                                                       2722, 5042, 1488, 1506,
                                                                       3112, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5477, 0, 3, 5042,
                                                                       2732, 5057, 1506, 1524,
                                                                       3142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5522, 0, 3, 5072,
                                                                       2812, 5117, 1560, 1596,
                                                                       3292, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5612, 0, 3, 5117,
                                                                       2842, 5162, 1596, 1632,
                                                                       3352, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5702, 0, 3, 5162,
                                                                       2872, 5207, 1632, 1668,
                                                                       3412, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5792, 0, 3, 5207,
                                                                       2902, 5252, 1668, 1704,
                                                                       3472, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5882, 0, 3, 5297,
                                                                       3022, 5342, 1776, 1812,
                                                                       3652, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5972, 0, 3, 5342,
                                                                       3052, 5387, 1812, 1848,
                                                                       3712, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6062, 0, 3, 5387,
                                                                       3082, 5432, 1848, 1884,
                                                                       3772, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6152, 0, 3, 5432,
                                                                       3112, 5477, 1884, 1920,
                                                                       3832, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 6242, 0, 3, 5522,
                                                                       3292, 5612, 1992, 2052,
                                                                       4092, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 6392, 0, 3, 5612,
                                                                       3352, 5702, 2052, 2112,
                                                                       4192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 6542, 0, 3, 5702,
                                                                       3412, 5792, 2112, 2172,
                                                                       4292, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 6692, 0, 3, 5882,
                                                                       3652, 5972, 2292, 2352,
                                                                       4592, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 6842, 0, 3, 5972,
                                                                       3712, 6062, 2352, 2412,
                                                                       4692, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 6992, 0, 3, 6062,
                                                                       3772, 6152, 2412, 2472,
                                                                       4792, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7142, 3, 2592,
                                                                       2602, 4892, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7163, 3, 2602,
                                                                       2612, 4907, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7184, 3, 2612,
                                                                       2622, 4922, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7205, 3, 2622,
                                                                       2632, 4937, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7226, 3, 2632,
                                                                       2642, 4952, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7247, 3, 2642,
                                                                       2652, 4967, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7268, 3, 2672,
                                                                       2682, 4982, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7289, 3, 2682,
                                                                       2692, 4997, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7310, 3, 2692,
                                                                       2702, 5012, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7331, 3, 2702,
                                                                       2712, 5027, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7352, 3, 2712,
                                                                       2722, 5042, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7373, 3, 2722,
                                                                       2732, 5057, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7394, 0, 3, 7142,
                                                                       4892, 7163, 2752, 2782,
                                                                       5072, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7457, 0, 3, 7163,
                                                                       4907, 7184, 2782, 2812,
                                                                       5117, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7520, 0, 3, 7184,
                                                                       4922, 7205, 2812, 2842,
                                                                       5162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7583, 0, 3, 7205,
                                                                       4937, 7226, 2842, 2872,
                                                                       5207, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7646, 0, 3, 7226,
                                                                       4952, 7247, 2872, 2902,
                                                                       5252, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7709, 0, 3, 7268,
                                                                       4982, 7289, 2962, 2992,
                                                                       5297, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7772, 0, 3, 7289,
                                                                       4997, 7310, 2992, 3022,
                                                                       5342, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7835, 0, 3, 7310,
                                                                       5012, 7331, 3022, 3052,
                                                                       5387, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7898, 0, 3, 7331,
                                                                       5027, 7352, 3052, 3082,
                                                                       5432, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7961, 0, 3, 7352,
                                                                       5042, 7373, 3082, 3112,
                                                                       5477, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 8024, 0, 3, 7394,
                                                                       5072, 7457, 3172, 3232,
                                                                       5522, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 8150, 0, 3, 7457,
                                                                       5117, 7520, 3232, 3292,
                                                                       5612, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 8276, 0, 3, 7520,
                                                                       5162, 7583, 3292, 3352,
                                                                       5702, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 8402, 0, 3, 7583,
                                                                       5207, 7646, 3352, 3412,
                                                                       5792, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 8528, 0, 3, 7709,
                                                                       5297, 7772, 3532, 3592,
                                                                       5882, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 8654, 0, 3, 7772,
                                                                       5342, 7835, 3592, 3652,
                                                                       5972, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 8780, 0, 3, 7835,
                                                                       5387, 7898, 3652, 3712,
                                                                       6062, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 8906, 0, 3, 7898,
                                                                       5432, 7961, 3712, 3772,
                                                                       6152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 9032, 0, 3, 8024,
                                                                       5522, 8150, 3892, 3992,
                                                                       6242, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 9242, 0, 3, 8150,
                                                                       5612, 8276, 3992, 4092,
                                                                       6392, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 9452, 0, 3, 8276,
                                                                       5702, 8402, 4092, 4192,
                                                                       6542, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 9662, 0, 3, 8528,
                                                                       5882, 8654, 4392, 4492,
                                                                       6692, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 9872, 0, 3, 8654,
                                                                       5972, 8780, 4492, 4592,
                                                                       6842, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 10082, 0, 3, 8780,
                                                                       6062, 8906, 4592, 4692,
                                                                       6992, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10292, 3, 4892,
                                                                       4907, 7184, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10320, 3, 4907,
                                                                       4922, 7205, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10348, 3, 4922,
                                                                       4937, 7226, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10376, 3, 4937,
                                                                       4952, 7247, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10404, 3, 4982,
                                                                       4997, 7310, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10432, 3, 4997,
                                                                       5012, 7331, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10460, 3, 5012,
                                                                       5027, 7352, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 10488, 3, 5027,
                                                                       5042, 7373, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 10516, 0, 3,
                                                                       10292, 7184, 10320, 5072,
                                                                       5117, 7520, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 10600, 0, 3,
                                                                       10320, 7205, 10348, 5117,
                                                                       5162, 7583, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 10684, 0, 3,
                                                                       10348, 7226, 10376, 5162,
                                                                       5207, 7646, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 10768, 0, 3,
                                                                       10404, 7310, 10432, 5297,
                                                                       5342, 7835, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 10852, 0, 3,
                                                                       10432, 7331, 10460, 5342,
                                                                       5387, 7898, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 10936, 0, 3,
                                                                       10460, 7352, 10488, 5387,
                                                                       5432, 7961, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 11020, 0, 3,
                                                                       10516, 7520, 10600, 5522,
                                                                       5612, 8276, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 11188, 0, 3,
                                                                       10600, 7583, 10684, 5612,
                                                                       5702, 8402, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 11356, 0, 3,
                                                                       10768, 7835, 10852, 5882,
                                                                       5972, 8780, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 11524, 0, 3,
                                                                       10852, 7898, 10936, 5972,
                                                                       6062, 8906, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 11692, 0, 3,
                                                                       11020, 8276, 11188, 6242,
                                                                       6392, 9452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 11972, 0, 3,
                                                                       11356, 8780, 11524, 6692,
                                                                       6842, 10082, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12252, 3, 7142,
                                                                       7163, 10292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12288, 3, 7163,
                                                                       7184, 10320, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12324, 3, 7184,
                                                                       7205, 10348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12360, 3, 7205,
                                                                       7226, 10376, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12396, 3, 7268,
                                                                       7289, 10404, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12432, 3, 7289,
                                                                       7310, 10432, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12468, 3, 7310,
                                                                       7331, 10460, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 12504, 3, 7331,
                                                                       7352, 10488, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 12540, 0, 3,
                                                                       12252, 10292, 12288, 7394,
                                                                       7457, 10516, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 12648, 0, 3,
                                                                       12288, 10320, 12324, 7457,
                                                                       7520, 10600, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 12756, 0, 3,
                                                                       12324, 10348, 12360, 7520,
                                                                       7583, 10684, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 12864, 0, 3,
                                                                       12396, 10404, 12432, 7709,
                                                                       7772, 10768, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 12972, 0, 3,
                                                                       12432, 10432, 12468, 7772,
                                                                       7835, 10852, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 13080, 0, 3,
                                                                       12468, 10460, 12504, 7835,
                                                                       7898, 10936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 13188, 0, 3,
                                                                       12540, 10516, 12648, 8024,
                                                                       8150, 11020, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 13404, 0, 3,
                                                                       12648, 10600, 12756, 8150,
                                                                       8276, 11188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 13620, 0, 3,
                                                                       12864, 10768, 12972, 8528,
                                                                       8654, 11356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 13836, 0, 3,
                                                                       12972, 10852, 13080, 8654,
                                                                       8780, 11524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 14052, 0, 3,
                                                                       13188, 11020, 13404, 9032,
                                                                       9242, 11692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 14412, 0, 3,
                                                                       13620, 11356, 13836, 9662,
                                                                       9872, 11972, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 14772, 13188, 216, ncols);

                    simdfunc::contract_primitives(buffer, 15078, 13620, 216, ncols);

                    simdfunc::contract_primitives(buffer, 15384, 14052, 360, ncols);

                    simdfunc::contract_primitives(buffer, 15894, 14412, 360, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 14988, 14772, 6, 1, nmax);

        simdtrf::transform_k_inner(buffer, 15294, 15078, 6, 1, nmax);

        simdtrf::transform_k_inner(buffer, 15744, 15384, 10, 1, nmax);

        simdtrf::transform_k_inner(buffer, 16254, 15894, 10, 1, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 16404, 14988, 15744, 15, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 16674, 15294, 16254, 15, nmax);

        simdtrf::transform_d_inner(buffer, 16944, 16674, 3, 15, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 16944, 75, nmax);

        simdtrf::transform_d_inner(buffer, 16944, 16404, 3, 15, nmax);

        simdtrf::transform_p_outer(values + 225 * nvalues + n * npairs, nvalues, buffer, 16944,
                                   75, nmax);
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
