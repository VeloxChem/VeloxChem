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


#include "SimdThreeCenterElectronRepulsionRecPGL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformL.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_pgl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_pgl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 38814, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 459 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 38814, 35358, 1875, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 13,
                                                             ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 61, 0, 3, 8, 9,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 67, 0, 3, 9, 10,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 73, 0, 3, 10, 11,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 79, 0, 3, 11, 12,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 85, 0, 3, 12, 13,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 91, 0, 3, 13, 14,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 97, 0, 3, 14, 15,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 103, 0, 3, 15, 16,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 109, 0, 3, 16, 17,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 115, 0, 3, 17, 18,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 121, 0, 3, 18, 19,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 127, 0, 3, 19, 20,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 133, 0, 3, 22, 25,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 143, 0, 3, 25, 28,
                                                                       67, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 28, 31,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 31, 34,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 34, 37,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 37, 40,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 40, 43,
                                                                       97, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 43, 46,
                                                                       103, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 46, 49,
                                                                       109, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 49, 52,
                                                                       115, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 233, 0, 3, 52, 55,
                                                                       121, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 61, 67,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 67, 73,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 73, 79,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 79, 85,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 303, 0, 3, 85, 91,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 91, 97,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 333, 0, 3, 97,
                                                                       103, 193, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 103,
                                                                       109, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 363, 0, 3, 109,
                                                                       115, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 115,
                                                                       121, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 393, 0, 3, 133,
                                                                       143, 243, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 414, 0, 3, 143,
                                                                       153, 258, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 435, 0, 3, 153,
                                                                       163, 273, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 456, 0, 3, 163,
                                                                       173, 288, 303, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 477, 0, 3, 173,
                                                                       183, 303, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 183,
                                                                       193, 318, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 519, 0, 3, 193,
                                                                       203, 333, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 540, 0, 3, 203,
                                                                       213, 348, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 561, 0, 3, 213,
                                                                       223, 363, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 582, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 585, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 588, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 591, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 594, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 597, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 600, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 603, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 606, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 609, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 612, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 615, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 618, 3, 10, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 627, 3, 11, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 636, 3, 12, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 645, 3, 13, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 654, 3, 14, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 663, 3, 15, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 672, 3, 16, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 681, 3, 17, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 690, 3, 18, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 699, 3, 19, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 708, 3, 20, 58,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 717, 3, 28, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 735, 3, 31, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 753, 3, 34, 85,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 771, 3, 37, 91,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 789, 3, 40, 97,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 807, 3, 43, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 825, 3, 46, 109,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 843, 3, 49, 115,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 861, 3, 52, 121,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 879, 3, 55, 127,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 897, 3, 73, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 927, 3, 79, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 957, 3, 85, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 987, 3, 91, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1017, 3, 97, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1047, 3, 103, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1077, 3, 109, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1107, 3, 115, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1137, 3, 121, 233,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1167, 3, 153, 273,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1212, 3, 163, 288,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1257, 3, 173, 303,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1302, 3, 183, 318,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1347, 3, 193, 333,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1392, 3, 203, 348,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1437, 3, 213, 363,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1482, 3, 223, 378,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1527, 3, 273, 435,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1590, 3, 288, 456,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1653, 3, 303, 477,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1716, 3, 318, 498,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1779, 3, 333, 519,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1842, 3, 348, 540,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1905, 3, 363, 561,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1968, 3, 8, 9,
                                                                       582, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1974, 3, 9, 10,
                                                                       585, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1980, 3, 10, 11,
                                                                       588, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1986, 3, 11, 12,
                                                                       591, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1992, 3, 12, 13,
                                                                       594, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1998, 3, 13, 14,
                                                                       597, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2004, 3, 14, 15,
                                                                       600, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2010, 3, 15, 16,
                                                                       603, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2016, 3, 16, 17,
                                                                       606, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2022, 3, 17, 18,
                                                                       609, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2028, 3, 18, 19,
                                                                       612, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2034, 3, 19, 20,
                                                                       615, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2040, 0, 3, 1968,
                                                                       582, 1974, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2058, 0, 3, 1974,
                                                                       585, 1980, 627, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2076, 0, 3, 1980,
                                                                       588, 1986, 636, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2094, 0, 3, 1986,
                                                                       591, 1992, 645, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2112, 0, 3, 1992,
                                                                       594, 1998, 654, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2130, 0, 3, 1998,
                                                                       597, 2004, 663, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2148, 0, 3, 2004,
                                                                       600, 2010, 672, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2166, 0, 3, 2010,
                                                                       603, 2016, 681, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2184, 0, 3, 2016,
                                                                       606, 2022, 690, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2202, 0, 3, 2022,
                                                                       609, 2028, 699, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2220, 0, 3, 2028,
                                                                       612, 2034, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2238, 0, 3, 2040,
                                                                       618, 2058, 61, 67, 717,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2274, 0, 3, 2058,
                                                                       627, 2076, 67, 73, 735,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2310, 0, 3, 2076,
                                                                       636, 2094, 73, 79, 753,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2346, 0, 3, 2094,
                                                                       645, 2112, 79, 85, 771,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2382, 0, 3, 2112,
                                                                       654, 2130, 85, 91, 789,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2418, 0, 3, 2130,
                                                                       663, 2148, 91, 97, 807,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2454, 0, 3, 2148,
                                                                       672, 2166, 97, 103, 825,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2490, 0, 3, 2166,
                                                                       681, 2184, 103, 109, 843,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2526, 0, 3, 2184,
                                                                       690, 2202, 109, 115, 861,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2562, 0, 3, 2202,
                                                                       699, 2220, 115, 121, 879,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2598, 0, 3, 2238,
                                                                       717, 2274, 133, 143, 897,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2658, 0, 3, 2274,
                                                                       735, 2310, 143, 153, 927,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2718, 0, 3, 2310,
                                                                       753, 2346, 153, 163, 957,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2778, 0, 3, 2346,
                                                                       771, 2382, 163, 173, 987,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2838, 0, 3, 2382,
                                                                       789, 2418, 173, 183, 1017,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2898, 0, 3, 2418,
                                                                       807, 2454, 183, 193, 1047,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2958, 0, 3, 2454,
                                                                       825, 2490, 193, 203, 1077,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3018, 0, 3, 2490,
                                                                       843, 2526, 203, 213, 1107,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3078, 0, 3, 2526,
                                                                       861, 2562, 213, 223, 1137,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3138, 0, 3, 2598,
                                                                       897, 2658, 243, 258, 1167,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3228, 0, 3, 2658,
                                                                       927, 2718, 258, 273, 1212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3318, 0, 3, 2718,
                                                                       957, 2778, 273, 288, 1257,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3408, 0, 3, 2778,
                                                                       987, 2838, 288, 303, 1302,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3498, 0, 3, 2838,
                                                                       1017, 2898, 303, 318,
                                                                       1347, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3588, 0, 3, 2898,
                                                                       1047, 2958, 318, 333,
                                                                       1392, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3678, 0, 3, 2958,
                                                                       1077, 3018, 333, 348,
                                                                       1437, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3768, 0, 3, 3018,
                                                                       1107, 3078, 348, 363,
                                                                       1482, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3858, 0, 3, 3138,
                                                                       1167, 3228, 393, 414,
                                                                       1527, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3984, 0, 3, 3228,
                                                                       1212, 3318, 414, 435,
                                                                       1590, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4110, 0, 3, 3318,
                                                                       1257, 3408, 435, 456,
                                                                       1653, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4236, 0, 3, 3408,
                                                                       1302, 3498, 456, 477,
                                                                       1716, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4362, 0, 3, 3498,
                                                                       1347, 3588, 477, 498,
                                                                       1779, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4488, 0, 3, 3588,
                                                                       1392, 3678, 498, 519,
                                                                       1842, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4614, 0, 3, 3678,
                                                                       1437, 3768, 519, 540,
                                                                       1905, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4740, 3, 582, 585,
                                                                       1980, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4750, 3, 585, 588,
                                                                       1986, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4760, 3, 588, 591,
                                                                       1992, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4770, 3, 591, 594,
                                                                       1998, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4780, 3, 594, 597,
                                                                       2004, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4790, 3, 597, 600,
                                                                       2010, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4800, 3, 600, 603,
                                                                       2016, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4810, 3, 603, 606,
                                                                       2022, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4820, 3, 606, 609,
                                                                       2028, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4830, 3, 609, 612,
                                                                       2034, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4840, 0, 3, 4740,
                                                                       1980, 4750, 2076, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4870, 0, 3, 4750,
                                                                       1986, 4760, 2094, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4900, 0, 3, 4760,
                                                                       1992, 4770, 2112, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4930, 0, 3, 4770,
                                                                       1998, 4780, 2130, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4960, 0, 3, 4780,
                                                                       2004, 4790, 2148, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4990, 0, 3, 4790,
                                                                       2010, 4800, 2166, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5020, 0, 3, 4800,
                                                                       2016, 4810, 2184, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5050, 0, 3, 4810,
                                                                       2022, 4820, 2202, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5080, 0, 3, 4820,
                                                                       2028, 4830, 2220, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5110, 0, 3, 4840,
                                                                       2076, 4870, 717, 735,
                                                                       2310, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5170, 0, 3, 4870,
                                                                       2094, 4900, 735, 753,
                                                                       2346, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5230, 0, 3, 4900,
                                                                       2112, 4930, 753, 771,
                                                                       2382, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5290, 0, 3, 4930,
                                                                       2130, 4960, 771, 789,
                                                                       2418, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5350, 0, 3, 4960,
                                                                       2148, 4990, 789, 807,
                                                                       2454, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5410, 0, 3, 4990,
                                                                       2166, 5020, 807, 825,
                                                                       2490, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5470, 0, 3, 5020,
                                                                       2184, 5050, 825, 843,
                                                                       2526, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5530, 0, 3, 5050,
                                                                       2202, 5080, 843, 861,
                                                                       2562, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5590, 0, 3, 5110,
                                                                       2310, 5170, 897, 927,
                                                                       2718, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5690, 0, 3, 5170,
                                                                       2346, 5230, 927, 957,
                                                                       2778, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5790, 0, 3, 5230,
                                                                       2382, 5290, 957, 987,
                                                                       2838, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5890, 0, 3, 5290,
                                                                       2418, 5350, 987, 1017,
                                                                       2898, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5990, 0, 3, 5350,
                                                                       2454, 5410, 1017, 1047,
                                                                       2958, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6090, 0, 3, 5410,
                                                                       2490, 5470, 1047, 1077,
                                                                       3018, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6190, 0, 3, 5470,
                                                                       2526, 5530, 1077, 1107,
                                                                       3078, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6290, 0, 3, 5590,
                                                                       2718, 5690, 1167, 1212,
                                                                       3318, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6440, 0, 3, 5690,
                                                                       2778, 5790, 1212, 1257,
                                                                       3408, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6590, 0, 3, 5790,
                                                                       2838, 5890, 1257, 1302,
                                                                       3498, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6740, 0, 3, 5890,
                                                                       2898, 5990, 1302, 1347,
                                                                       3588, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6890, 0, 3, 5990,
                                                                       2958, 6090, 1347, 1392,
                                                                       3678, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7040, 0, 3, 6090,
                                                                       3018, 6190, 1392, 1437,
                                                                       3768, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7190, 0, 3, 6290,
                                                                       3318, 6440, 1527, 1590,
                                                                       4110, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 6440,
                                                                       3408, 6590, 1590, 1653,
                                                                       4236, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7610, 0, 3, 6590,
                                                                       3498, 6740, 1653, 1716,
                                                                       4362, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7820, 0, 3, 6740,
                                                                       3588, 6890, 1716, 1779,
                                                                       4488, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8030, 0, 3, 6890,
                                                                       3678, 7040, 1779, 1842,
                                                                       4614, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8240, 3, 1968,
                                                                       1974, 4740, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8255, 3, 1974,
                                                                       1980, 4750, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8270, 3, 1980,
                                                                       1986, 4760, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8285, 3, 1986,
                                                                       1992, 4770, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8300, 3, 1992,
                                                                       1998, 4780, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8315, 3, 1998,
                                                                       2004, 4790, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8330, 3, 2004,
                                                                       2010, 4800, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8345, 3, 2010,
                                                                       2016, 4810, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8360, 3, 2016,
                                                                       2022, 4820, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8375, 3, 2022,
                                                                       2028, 4830, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8390, 0, 3, 8240,
                                                                       4740, 8255, 2040, 2058,
                                                                       4840, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8435, 0, 3, 8255,
                                                                       4750, 8270, 2058, 2076,
                                                                       4870, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8480, 0, 3, 8270,
                                                                       4760, 8285, 2076, 2094,
                                                                       4900, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8525, 0, 3, 8285,
                                                                       4770, 8300, 2094, 2112,
                                                                       4930, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8570, 0, 3, 8300,
                                                                       4780, 8315, 2112, 2130,
                                                                       4960, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8615, 0, 3, 8315,
                                                                       4790, 8330, 2130, 2148,
                                                                       4990, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8660, 0, 3, 8330,
                                                                       4800, 8345, 2148, 2166,
                                                                       5020, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8705, 0, 3, 8345,
                                                                       4810, 8360, 2166, 2184,
                                                                       5050, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8750, 0, 3, 8360,
                                                                       4820, 8375, 2184, 2202,
                                                                       5080, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8795, 0, 3, 8390,
                                                                       4840, 8435, 2238, 2274,
                                                                       5110, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8885, 0, 3, 8435,
                                                                       4870, 8480, 2274, 2310,
                                                                       5170, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8975, 0, 3, 8480,
                                                                       4900, 8525, 2310, 2346,
                                                                       5230, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9065, 0, 3, 8525,
                                                                       4930, 8570, 2346, 2382,
                                                                       5290, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9155, 0, 3, 8570,
                                                                       4960, 8615, 2382, 2418,
                                                                       5350, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9245, 0, 3, 8615,
                                                                       4990, 8660, 2418, 2454,
                                                                       5410, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9335, 0, 3, 8660,
                                                                       5020, 8705, 2454, 2490,
                                                                       5470, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9425, 0, 3, 8705,
                                                                       5050, 8750, 2490, 2526,
                                                                       5530, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9515, 0, 3, 8795,
                                                                       5110, 8885, 2598, 2658,
                                                                       5590, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9665, 0, 3, 8885,
                                                                       5170, 8975, 2658, 2718,
                                                                       5690, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9815, 0, 3, 8975,
                                                                       5230, 9065, 2718, 2778,
                                                                       5790, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9965, 0, 3, 9065,
                                                                       5290, 9155, 2778, 2838,
                                                                       5890, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10115, 0, 3, 9155,
                                                                       5350, 9245, 2838, 2898,
                                                                       5990, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10265, 0, 3, 9245,
                                                                       5410, 9335, 2898, 2958,
                                                                       6090, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10415, 0, 3, 9335,
                                                                       5470, 9425, 2958, 3018,
                                                                       6190, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 10565, 0, 3, 9515,
                                                                       5590, 9665, 3138, 3228,
                                                                       6290, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 10790, 0, 3, 9665,
                                                                       5690, 9815, 3228, 3318,
                                                                       6440, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 11015, 0, 3, 9815,
                                                                       5790, 9965, 3318, 3408,
                                                                       6590, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 11240, 0, 3, 9965,
                                                                       5890, 10115, 3408, 3498,
                                                                       6740, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 11465, 0, 3,
                                                                       10115, 5990, 10265, 3498,
                                                                       3588, 6890, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 11690, 0, 3,
                                                                       10265, 6090, 10415, 3588,
                                                                       3678, 7040, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 11915, 0, 3,
                                                                       10565, 6290, 10790, 3858,
                                                                       3984, 7190, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 12230, 0, 3,
                                                                       10790, 6440, 11015, 3984,
                                                                       4110, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 12545, 0, 3,
                                                                       11015, 6590, 11240, 4110,
                                                                       4236, 7610, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 12860, 0, 3,
                                                                       11240, 6740, 11465, 4236,
                                                                       4362, 7820, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 13175, 0, 3,
                                                                       11465, 6890, 11690, 4362,
                                                                       4488, 8030, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13490, 3, 4740,
                                                                       4750, 8270, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13511, 3, 4750,
                                                                       4760, 8285, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13532, 3, 4760,
                                                                       4770, 8300, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13553, 3, 4770,
                                                                       4780, 8315, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13574, 3, 4780,
                                                                       4790, 8330, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13595, 3, 4790,
                                                                       4800, 8345, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13616, 3, 4800,
                                                                       4810, 8360, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13637, 3, 4810,
                                                                       4820, 8375, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13658, 0, 3,
                                                                       13490, 8270, 13511, 4840,
                                                                       4870, 8480, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13721, 0, 3,
                                                                       13511, 8285, 13532, 4870,
                                                                       4900, 8525, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13784, 0, 3,
                                                                       13532, 8300, 13553, 4900,
                                                                       4930, 8570, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13847, 0, 3,
                                                                       13553, 8315, 13574, 4930,
                                                                       4960, 8615, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13910, 0, 3,
                                                                       13574, 8330, 13595, 4960,
                                                                       4990, 8660, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13973, 0, 3,
                                                                       13595, 8345, 13616, 4990,
                                                                       5020, 8705, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 14036, 0, 3,
                                                                       13616, 8360, 13637, 5020,
                                                                       5050, 8750, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14099, 0, 3,
                                                                       13658, 8480, 13721, 5110,
                                                                       5170, 8975, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14225, 0, 3,
                                                                       13721, 8525, 13784, 5170,
                                                                       5230, 9065, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14351, 0, 3,
                                                                       13784, 8570, 13847, 5230,
                                                                       5290, 9155, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14477, 0, 3,
                                                                       13847, 8615, 13910, 5290,
                                                                       5350, 9245, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14603, 0, 3,
                                                                       13910, 8660, 13973, 5350,
                                                                       5410, 9335, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14729, 0, 3,
                                                                       13973, 8705, 14036, 5410,
                                                                       5470, 9425, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 14855, 0, 3,
                                                                       14099, 8975, 14225, 5590,
                                                                       5690, 9815, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 15065, 0, 3,
                                                                       14225, 9065, 14351, 5690,
                                                                       5790, 9965, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 15275, 0, 3,
                                                                       14351, 9155, 14477, 5790,
                                                                       5890, 10115, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 15485, 0, 3,
                                                                       14477, 9245, 14603, 5890,
                                                                       5990, 10265, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 15695, 0, 3,
                                                                       14603, 9335, 14729, 5990,
                                                                       6090, 10415, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 15905, 0, 3,
                                                                       14855, 9815, 15065, 6290,
                                                                       6440, 11015, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 16220, 0, 3,
                                                                       15065, 9965, 15275, 6440,
                                                                       6590, 11240, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 16535, 0, 3,
                                                                       15275, 10115, 15485, 6590,
                                                                       6740, 11465, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 16850, 0, 3,
                                                                       15485, 10265, 15695, 6740,
                                                                       6890, 11690, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 17165, 0, 3,
                                                                       15905, 11015, 16220, 7190,
                                                                       7400, 12545, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 17606, 0, 3,
                                                                       16220, 11240, 16535, 7400,
                                                                       7610, 12860, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 18047, 0, 3,
                                                                       16535, 11465, 16850, 7610,
                                                                       7820, 13175, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18488, 3, 8240,
                                                                       8255, 13490, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18516, 3, 8255,
                                                                       8270, 13511, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18544, 3, 8270,
                                                                       8285, 13532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18572, 3, 8285,
                                                                       8300, 13553, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18600, 3, 8300,
                                                                       8315, 13574, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18628, 3, 8315,
                                                                       8330, 13595, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18656, 3, 8330,
                                                                       8345, 13616, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18684, 3, 8345,
                                                                       8360, 13637, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 18712, 0, 3,
                                                                       18488, 13490, 18516, 8390,
                                                                       8435, 13658, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 18796, 0, 3,
                                                                       18516, 13511, 18544, 8435,
                                                                       8480, 13721, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 18880, 0, 3,
                                                                       18544, 13532, 18572, 8480,
                                                                       8525, 13784, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 18964, 0, 3,
                                                                       18572, 13553, 18600, 8525,
                                                                       8570, 13847, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 19048, 0, 3,
                                                                       18600, 13574, 18628, 8570,
                                                                       8615, 13910, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 19132, 0, 3,
                                                                       18628, 13595, 18656, 8615,
                                                                       8660, 13973, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 19216, 0, 3,
                                                                       18656, 13616, 18684, 8660,
                                                                       8705, 14036, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 19300, 0, 3,
                                                                       18712, 13658, 18796, 8795,
                                                                       8885, 14099, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 19468, 0, 3,
                                                                       18796, 13721, 18880, 8885,
                                                                       8975, 14225, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 19636, 0, 3,
                                                                       18880, 13784, 18964, 8975,
                                                                       9065, 14351, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 19804, 0, 3,
                                                                       18964, 13847, 19048, 9065,
                                                                       9155, 14477, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 19972, 0, 3,
                                                                       19048, 13910, 19132, 9155,
                                                                       9245, 14603, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 20140, 0, 3,
                                                                       19132, 13973, 19216, 9245,
                                                                       9335, 14729, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 20308, 0, 3,
                                                                       19300, 14099, 19468, 9515,
                                                                       9665, 14855, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 20588, 0, 3,
                                                                       19468, 14225, 19636, 9665,
                                                                       9815, 15065, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 20868, 0, 3,
                                                                       19636, 14351, 19804, 9815,
                                                                       9965, 15275, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 21148, 0, 3,
                                                                       19804, 14477, 19972, 9965,
                                                                       10115, 15485, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 21428, 0, 3,
                                                                       19972, 14603, 20140,
                                                                       10115, 10265, 15695,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 21708, 0, 3,
                                                                       20308, 14855, 20588,
                                                                       10565, 10790, 15905,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 22128, 0, 3,
                                                                       20588, 15065, 20868,
                                                                       10790, 11015, 16220,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 22548, 0, 3,
                                                                       20868, 15275, 21148,
                                                                       11015, 11240, 16535,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 22968, 0, 3,
                                                                       21148, 15485, 21428,
                                                                       11240, 11465, 16850,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 23388, 0, 3,
                                                                       21708, 15905, 22128,
                                                                       11915, 12230, 17165,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 23976, 0, 3,
                                                                       22128, 16220, 22548,
                                                                       12230, 12545, 17606,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 24564, 0, 3,
                                                                       22548, 16535, 22968,
                                                                       12545, 12860, 18047,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25152, 3, 13490,
                                                                       13511, 18544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25188, 3, 13511,
                                                                       13532, 18572, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25224, 3, 13532,
                                                                       13553, 18600, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25260, 3, 13553,
                                                                       13574, 18628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25296, 3, 13574,
                                                                       13595, 18656, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25332, 3, 13595,
                                                                       13616, 18684, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 25368, 0, 3,
                                                                       25152, 18544, 25188,
                                                                       13658, 13721, 18880,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 25476, 0, 3,
                                                                       25188, 18572, 25224,
                                                                       13721, 13784, 18964,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 25584, 0, 3,
                                                                       25224, 18600, 25260,
                                                                       13784, 13847, 19048,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 25692, 0, 3,
                                                                       25260, 18628, 25296,
                                                                       13847, 13910, 19132,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 25800, 0, 3,
                                                                       25296, 18656, 25332,
                                                                       13910, 13973, 19216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 25908, 0, 3,
                                                                       25368, 18880, 25476,
                                                                       14099, 14225, 19636,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 26124, 0, 3,
                                                                       25476, 18964, 25584,
                                                                       14225, 14351, 19804,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 26340, 0, 3,
                                                                       25584, 19048, 25692,
                                                                       14351, 14477, 19972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 26556, 0, 3,
                                                                       25692, 19132, 25800,
                                                                       14477, 14603, 20140,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 26772, 0, 3,
                                                                       25908, 19636, 26124,
                                                                       14855, 15065, 20868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 27132, 0, 3,
                                                                       26124, 19804, 26340,
                                                                       15065, 15275, 21148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 27492, 0, 3,
                                                                       26340, 19972, 26556,
                                                                       15275, 15485, 21428,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 27852, 0, 3,
                                                                       26772, 20868, 27132,
                                                                       15905, 16220, 22548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 28392, 0, 3,
                                                                       27132, 21148, 27492,
                                                                       16220, 16535, 22968,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 28932, 0, 3,
                                                                       27852, 22548, 28392,
                                                                       17165, 17606, 24564,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29688, 3, 18488,
                                                                       18516, 25152, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29733, 3, 18516,
                                                                       18544, 25188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29778, 3, 18544,
                                                                       18572, 25224, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29823, 3, 18572,
                                                                       18600, 25260, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29868, 3, 18600,
                                                                       18628, 25296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29913, 3, 18628,
                                                                       18656, 25332, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 29958, 0, 3,
                                                                       29688, 25152, 29733,
                                                                       18712, 18796, 25368,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 30093, 0, 3,
                                                                       29733, 25188, 29778,
                                                                       18796, 18880, 25476,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 30228, 0, 3,
                                                                       29778, 25224, 29823,
                                                                       18880, 18964, 25584,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 30363, 0, 3,
                                                                       29823, 25260, 29868,
                                                                       18964, 19048, 25692,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 30498, 0, 3,
                                                                       29868, 25296, 29913,
                                                                       19048, 19132, 25800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 30633, 0, 3,
                                                                       29958, 25368, 30093,
                                                                       19300, 19468, 25908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 30903, 0, 3,
                                                                       30093, 25476, 30228,
                                                                       19468, 19636, 26124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 31173, 0, 3,
                                                                       30228, 25584, 30363,
                                                                       19636, 19804, 26340,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 31443, 0, 3,
                                                                       30363, 25692, 30498,
                                                                       19804, 19972, 26556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 31713, 0, 3,
                                                                       30633, 25908, 30903,
                                                                       20308, 20588, 26772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 32163, 0, 3,
                                                                       30903, 26124, 31173,
                                                                       20588, 20868, 27132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 32613, 0, 3,
                                                                       31173, 26340, 31443,
                                                                       20868, 21148, 27492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 33063, 0, 3,
                                                                       31713, 26772, 32163,
                                                                       21708, 22128, 27852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 33738, 0, 3,
                                                                       32163, 27132, 32613,
                                                                       22128, 22548, 28392,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 34413, 0, 3,
                                                                       33063, 27852, 33738,
                                                                       23388, 23976, 28932,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 35358, 33063, 675, ncols);

                    simdfunc::contract_primitives(buffer, 36288, 34413, 945, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 36033, 35358, 15, 1, nmax);

        simdtrf::transform_l_inner(buffer, 37233, 36288, 21, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 37590, 36033, 37233, 17, nmax);

        simdtrf::transform_g_inner(buffer, 38355, 37590, 3, 17, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 38355, 153, nmax);
    }

    for (size_t m = 0; m < 459; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
