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


#include "SimdThreeCenterElectronRepulsionRecPDK.hpp"

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
compute_pdk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_pdk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 8701, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 225 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 8701, 7390, 666, dimensions);

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
                                                        5, 6, 7, 8, 9, 10}, ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 18, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 45, 0, 3, 8, 9,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 51, 0, 3, 9, 10,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 57, 0, 3, 10, 11,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 63, 0, 3, 11, 12,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 69, 0, 3, 12, 13,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 75, 0, 3, 13, 14,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 81, 0, 3, 14, 15,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 87, 0, 3, 15, 16,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 93, 0, 3, 18, 21,
                                                                       45, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 103, 0, 3, 21, 24,
                                                                       51, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 113, 0, 3, 24, 27,
                                                                       57, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 123, 0, 3, 27, 30,
                                                                       63, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 133, 0, 3, 30, 33,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 143, 0, 3, 33, 36,
                                                                       75, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 36, 39,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 163, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 166, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 169, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 172, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 175, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 178, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 181, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 184, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 187, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 190, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 193, 3, 10, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 202, 3, 11, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 211, 3, 12, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 220, 3, 13, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 229, 3, 14, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 238, 3, 15, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 247, 3, 16, 42,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 256, 3, 18, 45,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 274, 3, 21, 51,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 292, 3, 24, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 310, 3, 27, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 328, 3, 30, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 346, 3, 33, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 364, 3, 36, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 382, 3, 39, 87,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 400, 3, 45, 93,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 430, 3, 51, 103,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 460, 3, 57, 113,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 490, 3, 63, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 520, 3, 69, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 550, 3, 75, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 580, 3, 81, 153,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 610, 3, 8, 9, 169,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 616, 3, 9, 10,
                                                                       172, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 622, 3, 10, 11,
                                                                       175, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 628, 3, 11, 12,
                                                                       178, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 634, 3, 12, 13,
                                                                       181, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 640, 3, 13, 14,
                                                                       184, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 646, 3, 14, 15,
                                                                       187, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 652, 3, 15, 16,
                                                                       190, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 658, 0, 3, 610,
                                                                       169, 616, 193, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 676, 0, 3, 616,
                                                                       172, 622, 202, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 694, 0, 3, 622,
                                                                       175, 628, 211, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 712, 0, 3, 628,
                                                                       178, 634, 220, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 730, 0, 3, 634,
                                                                       181, 640, 229, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 748, 0, 3, 640,
                                                                       184, 646, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 766, 0, 3, 646,
                                                                       187, 652, 247, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 784, 0, 3, 658,
                                                                       193, 676, 45, 51, 292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 820, 0, 3, 676,
                                                                       202, 694, 51, 57, 310,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 856, 0, 3, 694,
                                                                       211, 712, 57, 63, 328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 892, 0, 3, 712,
                                                                       220, 730, 63, 69, 346,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 928, 0, 3, 730,
                                                                       229, 748, 69, 75, 364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 964, 0, 3, 748,
                                                                       238, 766, 75, 81, 382,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1000, 0, 3, 784,
                                                                       292, 820, 93, 103, 460,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1060, 0, 3, 820,
                                                                       310, 856, 103, 113, 490,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1120, 0, 3, 856,
                                                                       328, 892, 113, 123, 520,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1180, 0, 3, 892,
                                                                       346, 928, 123, 133, 550,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 928,
                                                                       364, 964, 133, 143, 580,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1300, 3, 163, 166,
                                                                       610, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1310, 3, 166, 169,
                                                                       616, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1320, 3, 169, 172,
                                                                       622, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1330, 3, 172, 175,
                                                                       628, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1340, 3, 175, 178,
                                                                       634, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1350, 3, 178, 181,
                                                                       640, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1360, 3, 181, 184,
                                                                       646, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1370, 3, 184, 187,
                                                                       652, ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 1300,
                                                                       610, 1310, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1410, 0, 3, 1310,
                                                                       616, 1320, 676, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1440, 0, 3, 1320,
                                                                       622, 1330, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1470, 0, 3, 1330,
                                                                       628, 1340, 712, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1500, 0, 3, 1340,
                                                                       634, 1350, 730, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1530, 0, 3, 1350,
                                                                       640, 1360, 748, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1560, 0, 3, 1360,
                                                                       646, 1370, 766, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1590, 0, 3, 1380,
                                                                       658, 1410, 256, 274, 784,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1650, 0, 3, 1410,
                                                                       676, 1440, 274, 292, 820,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1710, 0, 3, 1440,
                                                                       694, 1470, 292, 310, 856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1770, 0, 3, 1470,
                                                                       712, 1500, 310, 328, 892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1830, 0, 3, 1500,
                                                                       730, 1530, 328, 346, 928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1890, 0, 3, 1530,
                                                                       748, 1560, 346, 364, 964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 1950, 0, 3, 1590,
                                                                       784, 1650, 400, 430, 1000,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 2050, 0, 3, 1650,
                                                                       820, 1710, 430, 460, 1060,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 2150, 0, 3, 1710,
                                                                       856, 1770, 460, 490, 1120,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 2250, 0, 3, 1770,
                                                                       892, 1830, 490, 520, 1180,
                                                                       ncols, gamma, p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 2350, 0, 3, 1830,
                                                                       928, 1890, 520, 550, 1240,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2450, 3, 610, 616,
                                                                       1320, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2465, 3, 616, 622,
                                                                       1330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2480, 3, 622, 628,
                                                                       1340, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2495, 3, 628, 634,
                                                                       1350, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2510, 3, 634, 640,
                                                                       1360, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2525, 3, 640, 646,
                                                                       1370, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2540, 0, 3, 2450,
                                                                       1320, 2465, 658, 676,
                                                                       1440, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2585, 0, 3, 2465,
                                                                       1330, 2480, 676, 694,
                                                                       1470, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2630, 0, 3, 2480,
                                                                       1340, 2495, 694, 712,
                                                                       1500, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2675, 0, 3, 2495,
                                                                       1350, 2510, 712, 730,
                                                                       1530, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2720, 0, 3, 2510,
                                                                       1360, 2525, 730, 748,
                                                                       1560, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 2765, 0, 3, 2540,
                                                                       1440, 2585, 784, 820,
                                                                       1710, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 2855, 0, 3, 2585,
                                                                       1470, 2630, 820, 856,
                                                                       1770, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 2945, 0, 3, 2630,
                                                                       1500, 2675, 856, 892,
                                                                       1830, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 3035, 0, 3, 2675,
                                                                       1530, 2720, 892, 928,
                                                                       1890, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 3125, 0, 3, 2765,
                                                                       1710, 2855, 1000, 1060,
                                                                       2150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 3275, 0, 3, 2855,
                                                                       1770, 2945, 1060, 1120,
                                                                       2250, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 3425, 0, 3, 2945,
                                                                       1830, 3035, 1120, 1180,
                                                                       2350, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3575, 3, 1300,
                                                                       1310, 2450, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3596, 3, 1310,
                                                                       1320, 2465, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3617, 3, 1320,
                                                                       1330, 2480, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3638, 3, 1330,
                                                                       1340, 2495, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3659, 3, 1340,
                                                                       1350, 2510, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3680, 3, 1350,
                                                                       1360, 2525, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 3701, 0, 3, 3575,
                                                                       2450, 3596, 1380, 1410,
                                                                       2540, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 3764, 0, 3, 3596,
                                                                       2465, 3617, 1410, 1440,
                                                                       2585, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 3827, 0, 3, 3617,
                                                                       2480, 3638, 1440, 1470,
                                                                       2630, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 3890, 0, 3, 3638,
                                                                       2495, 3659, 1470, 1500,
                                                                       2675, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 3953, 0, 3, 3659,
                                                                       2510, 3680, 1500, 1530,
                                                                       2720, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 4016, 0, 3, 3701,
                                                                       2540, 3764, 1590, 1650,
                                                                       2765, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 4142, 0, 3, 3764,
                                                                       2585, 3827, 1650, 1710,
                                                                       2855, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 4268, 0, 3, 3827,
                                                                       2630, 3890, 1710, 1770,
                                                                       2945, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 4394, 0, 3, 3890,
                                                                       2675, 3953, 1770, 1830,
                                                                       3035, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 4520, 0, 3, 4016,
                                                                       2765, 4142, 1950, 2050,
                                                                       3125, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 4730, 0, 3, 4142,
                                                                       2855, 4268, 2050, 2150,
                                                                       3275, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 4940, 0, 3, 4268,
                                                                       2945, 4394, 2150, 2250,
                                                                       3425, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5150, 3, 2450,
                                                                       2465, 3617, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5178, 3, 2465,
                                                                       2480, 3638, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5206, 3, 2480,
                                                                       2495, 3659, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5234, 3, 2495,
                                                                       2510, 3680, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 5262, 0, 3, 5150,
                                                                       3617, 5178, 2540, 2585,
                                                                       3827, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 5346, 0, 3, 5178,
                                                                       3638, 5206, 2585, 2630,
                                                                       3890, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 5430, 0, 3, 5206,
                                                                       3659, 5234, 2630, 2675,
                                                                       3953, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 5514, 0, 3, 5262,
                                                                       3827, 5346, 2765, 2855,
                                                                       4268, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 5682, 0, 3, 5346,
                                                                       3890, 5430, 2855, 2945,
                                                                       4394, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 5850, 0, 3, 5514,
                                                                       4268, 5682, 3125, 3275,
                                                                       4940, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 6130, 3, 3575,
                                                                       3596, 5150, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 6166, 3, 3596,
                                                                       3617, 5178, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 6202, 3, 3617,
                                                                       3638, 5206, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 6238, 3, 3638,
                                                                       3659, 5234, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 6274, 0, 3, 6130,
                                                                       5150, 6166, 3701, 3764,
                                                                       5262, ncols, gamma, p,
                                                                       q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 6382, 0, 3, 6166,
                                                                       5178, 6202, 3764, 3827,
                                                                       5346, ncols, gamma, p,
                                                                       q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 6490, 0, 3, 6202,
                                                                       5206, 6238, 3827, 3890,
                                                                       5430, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 6598, 0, 3, 6274,
                                                                       5262, 6382, 4016, 4142,
                                                                       5514, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 6814, 0, 3, 6382,
                                                                       5346, 6490, 4142, 4268,
                                                                       5682, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 7030, 0, 3, 6598,
                                                                       5514, 6814, 4520, 4730,
                                                                       5850, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 7390, 6598, 216, ncols);

                    simdfunc::contract_primitives(buffer, 7696, 7030, 360, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 7606, 7390, 6, 1, nmax);

        simdtrf::transform_k_inner(buffer, 8056, 7696, 10, 1, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 8206, 7606, 8056, 15, nmax);

        simdtrf::transform_d_inner(buffer, 8476, 8206, 3, 15, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 8476, 75, nmax);
    }

    for (size_t m = 0; m < 225; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
