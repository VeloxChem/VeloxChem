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


#include "SimdThreeCenterElectronRepulsionRecSFL.hpp"

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
#include "SimdTransformF.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_sfl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_sfl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 11603, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 119 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 11603, 10983, 450, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 11,
                                                             ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 53, 0, 3, 8, 9,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 59, 0, 3, 9, 10,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 65, 0, 3, 10, 11,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 71, 0, 3, 11, 12,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 77, 0, 3, 12, 13,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 83, 0, 3, 13, 14,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 89, 0, 3, 14, 15,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 95, 0, 3, 15, 16,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 101, 0, 3, 16, 17,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 107, 0, 3, 17, 18,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 113, 0, 3, 20, 23,
                                                                       53, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 123, 0, 3, 23, 26,
                                                                       59, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 133, 0, 3, 26, 29,
                                                                       65, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 143, 0, 3, 29, 32,
                                                                       71, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 32, 35,
                                                                       77, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 35, 38,
                                                                       83, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 38, 41,
                                                                       89, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 41, 44,
                                                                       95, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 44, 47,
                                                                       101, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 203, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 206, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 209, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 212, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 215, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 218, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 221, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 224, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 227, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 230, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 233, 3, 10, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 242, 3, 11, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 251, 3, 12, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 260, 3, 13, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 269, 3, 14, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 278, 3, 15, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 287, 3, 16, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 296, 3, 17, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 305, 3, 18, 50,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 314, 3, 26, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 332, 3, 29, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 350, 3, 32, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 368, 3, 35, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 386, 3, 38, 89,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 404, 3, 41, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 422, 3, 44, 101,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 440, 3, 47, 107,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 458, 3, 65, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 488, 3, 71, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 518, 3, 77, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 548, 3, 83, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 578, 3, 89, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 608, 3, 95, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 638, 3, 101, 193,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 668, 3, 8, 9, 203,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 674, 3, 9, 10,
                                                                       206, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 680, 3, 10, 11,
                                                                       209, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 686, 3, 11, 12,
                                                                       212, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 692, 3, 12, 13,
                                                                       215, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 698, 3, 13, 14,
                                                                       218, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 704, 3, 14, 15,
                                                                       221, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 710, 3, 15, 16,
                                                                       224, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 716, 3, 16, 17,
                                                                       227, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 722, 3, 17, 18,
                                                                       230, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 728, 0, 3, 668,
                                                                       203, 674, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 746, 0, 3, 674,
                                                                       206, 680, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 764, 0, 3, 680,
                                                                       209, 686, 251, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 782, 0, 3, 686,
                                                                       212, 692, 260, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 800, 0, 3, 692,
                                                                       215, 698, 269, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 818, 0, 3, 698,
                                                                       218, 704, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 836, 0, 3, 704,
                                                                       221, 710, 287, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 854, 0, 3, 710,
                                                                       224, 716, 296, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 872, 0, 3, 716,
                                                                       227, 722, 305, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 890, 0, 3, 728,
                                                                       233, 746, 53, 59, 314,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 926, 0, 3, 746,
                                                                       242, 764, 59, 65, 332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 962, 0, 3, 764,
                                                                       251, 782, 65, 71, 350,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 998, 0, 3, 782,
                                                                       260, 800, 71, 77, 368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1034, 0, 3, 800,
                                                                       269, 818, 77, 83, 386,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1070, 0, 3, 818,
                                                                       278, 836, 83, 89, 404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1106, 0, 3, 836,
                                                                       287, 854, 89, 95, 422,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 854,
                                                                       296, 872, 95, 101, 440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1178, 0, 3, 890,
                                                                       314, 926, 113, 123, 458,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1238, 0, 3, 926,
                                                                       332, 962, 123, 133, 488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1298, 0, 3, 962,
                                                                       350, 998, 133, 143, 518,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1358, 0, 3, 998,
                                                                       368, 1034, 143, 153, 548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1418, 0, 3, 1034,
                                                                       386, 1070, 153, 163, 578,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1478, 0, 3, 1070,
                                                                       404, 1106, 163, 173, 608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1538, 0, 3, 1106,
                                                                       422, 1142, 173, 183, 638,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1598, 3, 203, 206,
                                                                       680, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1608, 3, 206, 209,
                                                                       686, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1618, 3, 209, 212,
                                                                       692, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1628, 3, 212, 215,
                                                                       698, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1638, 3, 215, 218,
                                                                       704, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1648, 3, 218, 221,
                                                                       710, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1658, 3, 221, 224,
                                                                       716, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1668, 3, 224, 227,
                                                                       722, ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1678, 0, 3, 1598,
                                                                       680, 1608, 764, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1708, 0, 3, 1608,
                                                                       686, 1618, 782, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1738, 0, 3, 1618,
                                                                       692, 1628, 800, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1768, 0, 3, 1628,
                                                                       698, 1638, 818, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1798, 0, 3, 1638,
                                                                       704, 1648, 836, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 1648,
                                                                       710, 1658, 854, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1858, 0, 3, 1658,
                                                                       716, 1668, 872, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1888, 0, 3, 1678,
                                                                       764, 1708, 314, 332, 962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1948, 0, 3, 1708,
                                                                       782, 1738, 332, 350, 998,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2008, 0, 3, 1738,
                                                                       800, 1768, 350, 368, 1034,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2068, 0, 3, 1768,
                                                                       818, 1798, 368, 386, 1070,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2128, 0, 3, 1798,
                                                                       836, 1828, 386, 404, 1106,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2188, 0, 3, 1828,
                                                                       854, 1858, 404, 422, 1142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 2248, 0, 3, 1888,
                                                                       962, 1948, 458, 488, 1298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 2348, 0, 3, 1948,
                                                                       998, 2008, 488, 518, 1358,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 2448, 0, 3, 2008,
                                                                       1034, 2068, 518, 548,
                                                                       1418, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 2548, 0, 3, 2068,
                                                                       1070, 2128, 548, 578,
                                                                       1478, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 2128,
                                                                       1106, 2188, 578, 608,
                                                                       1538, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2748, 3, 668, 674,
                                                                       1598, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2763, 3, 674, 680,
                                                                       1608, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2778, 3, 680, 686,
                                                                       1618, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2793, 3, 686, 692,
                                                                       1628, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2808, 3, 692, 698,
                                                                       1638, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2823, 3, 698, 704,
                                                                       1648, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2838, 3, 704, 710,
                                                                       1658, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2853, 3, 710, 716,
                                                                       1668, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2868, 0, 3, 2748,
                                                                       1598, 2763, 728, 746,
                                                                       1678, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2913, 0, 3, 2763,
                                                                       1608, 2778, 746, 764,
                                                                       1708, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2958, 0, 3, 2778,
                                                                       1618, 2793, 764, 782,
                                                                       1738, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3003, 0, 3, 2793,
                                                                       1628, 2808, 782, 800,
                                                                       1768, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3048, 0, 3, 2808,
                                                                       1638, 2823, 800, 818,
                                                                       1798, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3093, 0, 3, 2823,
                                                                       1648, 2838, 818, 836,
                                                                       1828, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3138, 0, 3, 2838,
                                                                       1658, 2853, 836, 854,
                                                                       1858, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 3183, 0, 3, 2868,
                                                                       1678, 2913, 890, 926,
                                                                       1888, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 3273, 0, 3, 2913,
                                                                       1708, 2958, 926, 962,
                                                                       1948, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 3363, 0, 3, 2958,
                                                                       1738, 3003, 962, 998,
                                                                       2008, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 3453, 0, 3, 3003,
                                                                       1768, 3048, 998, 1034,
                                                                       2068, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 3543, 0, 3, 3048,
                                                                       1798, 3093, 1034, 1070,
                                                                       2128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 3633, 0, 3, 3093,
                                                                       1828, 3138, 1070, 1106,
                                                                       2188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 3723, 0, 3, 3183,
                                                                       1888, 3273, 1178, 1238,
                                                                       2248, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 3873, 0, 3, 3273,
                                                                       1948, 3363, 1238, 1298,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 4023, 0, 3, 3363,
                                                                       2008, 3453, 1298, 1358,
                                                                       2448, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 4173, 0, 3, 3453,
                                                                       2068, 3543, 1358, 1418,
                                                                       2548, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 4323, 0, 3, 3543,
                                                                       2128, 3633, 1418, 1478,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4473, 3, 1598,
                                                                       1608, 2778, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4494, 3, 1608,
                                                                       1618, 2793, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4515, 3, 1618,
                                                                       1628, 2808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4536, 3, 1628,
                                                                       1638, 2823, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4557, 3, 1638,
                                                                       1648, 2838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4578, 3, 1648,
                                                                       1658, 2853, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 4599, 0, 3, 4473,
                                                                       2778, 4494, 1678, 1708,
                                                                       2958, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 4662, 0, 3, 4494,
                                                                       2793, 4515, 1708, 1738,
                                                                       3003, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 4725, 0, 3, 4515,
                                                                       2808, 4536, 1738, 1768,
                                                                       3048, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 4788, 0, 3, 4536,
                                                                       2823, 4557, 1768, 1798,
                                                                       3093, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 4851, 0, 3, 4557,
                                                                       2838, 4578, 1798, 1828,
                                                                       3138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 4914, 0, 3, 4599,
                                                                       2958, 4662, 1888, 1948,
                                                                       3363, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 5040, 0, 3, 4662,
                                                                       3003, 4725, 1948, 2008,
                                                                       3453, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 5166, 0, 3, 4725,
                                                                       3048, 4788, 2008, 2068,
                                                                       3543, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 5292, 0, 3, 4788,
                                                                       3093, 4851, 2068, 2128,
                                                                       3633, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 5418, 0, 3, 4914,
                                                                       3363, 5040, 2248, 2348,
                                                                       4023, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 5628, 0, 3, 5040,
                                                                       3453, 5166, 2348, 2448,
                                                                       4173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 5838, 0, 3, 5166,
                                                                       3543, 5292, 2448, 2548,
                                                                       4323, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 6048, 3, 2748,
                                                                       2763, 4473, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 6076, 3, 2763,
                                                                       2778, 4494, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 6104, 3, 2778,
                                                                       2793, 4515, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 6132, 3, 2793,
                                                                       2808, 4536, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 6160, 3, 2808,
                                                                       2823, 4557, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 6188, 3, 2823,
                                                                       2838, 4578, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 6216, 0, 3, 6048,
                                                                       4473, 6076, 2868, 2913,
                                                                       4599, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 6300, 0, 3, 6076,
                                                                       4494, 6104, 2913, 2958,
                                                                       4662, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 6384, 0, 3, 6104,
                                                                       4515, 6132, 2958, 3003,
                                                                       4725, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 6468, 0, 3, 6132,
                                                                       4536, 6160, 3003, 3048,
                                                                       4788, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 6552, 0, 3, 6160,
                                                                       4557, 6188, 3048, 3093,
                                                                       4851, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 6636, 0, 3, 6216,
                                                                       4599, 6300, 3183, 3273,
                                                                       4914, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 6804, 0, 3, 6300,
                                                                       4662, 6384, 3273, 3363,
                                                                       5040, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 6972, 0, 3, 6384,
                                                                       4725, 6468, 3363, 3453,
                                                                       5166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 7140, 0, 3, 6468,
                                                                       4788, 6552, 3453, 3543,
                                                                       5292, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 7308, 0, 3, 6636,
                                                                       4914, 6804, 3723, 3873,
                                                                       5418, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 7588, 0, 3, 6804,
                                                                       5040, 6972, 3873, 4023,
                                                                       5628, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 7868, 0, 3, 6972,
                                                                       5166, 7140, 4023, 4173,
                                                                       5838, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 8148, 3, 4473,
                                                                       4494, 6104, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 8184, 3, 4494,
                                                                       4515, 6132, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 8220, 3, 4515,
                                                                       4536, 6160, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 8256, 3, 4536,
                                                                       4557, 6188, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 8292, 0, 3, 8148,
                                                                       6104, 8184, 4599, 4662,
                                                                       6384, ncols, gamma, p,
                                                                       q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 8400, 0, 3, 8184,
                                                                       6132, 8220, 4662, 4725,
                                                                       6468, ncols, gamma, p,
                                                                       q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 8508, 0, 3, 8220,
                                                                       6160, 8256, 4725, 4788,
                                                                       6552, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 8616, 0, 3, 8292,
                                                                       6384, 8400, 4914, 5040,
                                                                       6972, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 8832, 0, 3, 8400,
                                                                       6468, 8508, 5040, 5166,
                                                                       7140, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 9048, 0, 3, 8616,
                                                                       6972, 8832, 5418, 5628,
                                                                       7868, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 9408, 3, 6048,
                                                                       6076, 8148, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 9453, 3, 6076,
                                                                       6104, 8184, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 9498, 3, 6104,
                                                                       6132, 8220, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 9543, 3, 6132,
                                                                       6160, 8256, ncols, gamma,
                                                                       p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 9588, 0, 3, 9408,
                                                                       8148, 9453, 6216, 6300,
                                                                       8292, ncols, gamma, p,
                                                                       q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 9723, 0, 3, 9453,
                                                                       8184, 9498, 6300, 6384,
                                                                       8400, ncols, gamma, p,
                                                                       q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 9858, 0, 3, 9498,
                                                                       8220, 9543, 6384, 6468,
                                                                       8508, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 9993, 0, 3, 9588,
                                                                       8292, 9723, 6636, 6804,
                                                                       8616, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 10263, 0, 3, 9723,
                                                                       8400, 9858, 6804, 6972,
                                                                       8832, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 10533, 0, 3, 9993,
                                                                       8616, 10263, 7308, 7588,
                                                                       9048, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 10983, 10533, 450, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 11433, 10983, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 11433, 17, nmax);
    }

    for (size_t m = 0; m < 119; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
