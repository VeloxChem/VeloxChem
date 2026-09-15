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


#include "SimdThreeCenterElectronRepulsionRecPGG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_pgg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_pgg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 7162, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 243 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 7162, 5650, 675, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 9, ncols,
                                                             fj, 6, fq);

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

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 45, 51,
                                                                       93, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 178, 0, 3, 51, 57,
                                                                       103, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 57, 63,
                                                                       113, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 208, 0, 3, 63, 69,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 69, 75,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 75, 81,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 253, 0, 3, 93,
                                                                       103, 163, 178, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 274, 0, 3, 103,
                                                                       113, 178, 193, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 295, 0, 3, 113,
                                                                       123, 193, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 316, 0, 3, 123,
                                                                       133, 208, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 337, 0, 3, 133,
                                                                       143, 223, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 358, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 361, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 364, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 367, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 370, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 373, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 376, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 379, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 382, 3, 10, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 391, 3, 11, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 400, 3, 12, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 409, 3, 13, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 418, 3, 14, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 427, 3, 15, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 436, 3, 16, 42,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 445, 3, 24, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 463, 3, 27, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 481, 3, 30, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 499, 3, 33, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 517, 3, 36, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 535, 3, 39, 87,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 553, 3, 57, 113,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 583, 3, 63, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 613, 3, 69, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 643, 3, 75, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 673, 3, 81, 153,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 703, 3, 113, 193,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 748, 3, 123, 208,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 793, 3, 133, 223,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 838, 3, 143, 238,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 883, 3, 193, 295,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 946, 3, 208, 316,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1009, 3, 223, 337,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1072, 3, 8, 9,
                                                                       358, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1078, 3, 9, 10,
                                                                       361, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1084, 3, 10, 11,
                                                                       364, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1090, 3, 11, 12,
                                                                       367, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1096, 3, 12, 13,
                                                                       370, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1102, 3, 13, 14,
                                                                       373, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1108, 3, 14, 15,
                                                                       376, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1114, 3, 15, 16,
                                                                       379, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1120, 0, 3, 1072,
                                                                       358, 1078, 382, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1138, 0, 3, 1078,
                                                                       361, 1084, 391, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 1084,
                                                                       364, 1090, 400, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1174, 0, 3, 1090,
                                                                       367, 1096, 409, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1192, 0, 3, 1096,
                                                                       370, 1102, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1210, 0, 3, 1102,
                                                                       373, 1108, 427, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1228, 0, 3, 1108,
                                                                       376, 1114, 436, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1246, 0, 3, 1120,
                                                                       382, 1138, 45, 51, 445,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1282, 0, 3, 1138,
                                                                       391, 1156, 51, 57, 463,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1318, 0, 3, 1156,
                                                                       400, 1174, 57, 63, 481,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1354, 0, 3, 1174,
                                                                       409, 1192, 63, 69, 499,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1390, 0, 3, 1192,
                                                                       418, 1210, 69, 75, 517,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1426, 0, 3, 1210,
                                                                       427, 1228, 75, 81, 535,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1462, 0, 3, 1246,
                                                                       445, 1282, 93, 103, 553,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1522, 0, 3, 1282,
                                                                       463, 1318, 103, 113, 583,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1582, 0, 3, 1318,
                                                                       481, 1354, 113, 123, 613,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1642, 0, 3, 1354,
                                                                       499, 1390, 123, 133, 643,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1702, 0, 3, 1390,
                                                                       517, 1426, 133, 143, 673,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 1762, 0, 3, 1462,
                                                                       553, 1522, 163, 178, 703,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 1852, 0, 3, 1522,
                                                                       583, 1582, 178, 193, 748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 1942, 0, 3, 1582,
                                                                       613, 1642, 193, 208, 793,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2032, 0, 3, 1642,
                                                                       643, 1702, 208, 223, 838,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2122, 0, 3, 1762,
                                                                       703, 1852, 253, 274, 883,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2248, 0, 3, 1852,
                                                                       748, 1942, 274, 295, 946,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2374, 0, 3, 1942,
                                                                       793, 2032, 295, 316, 1009,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2500, 3, 358, 361,
                                                                       1084, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2510, 3, 361, 364,
                                                                       1090, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2520, 3, 364, 367,
                                                                       1096, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2530, 3, 367, 370,
                                                                       1102, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2540, 3, 370, 373,
                                                                       1108, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2550, 3, 373, 376,
                                                                       1114, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2560, 0, 3, 2500,
                                                                       1084, 2510, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2590, 0, 3, 2510,
                                                                       1090, 2520, 1174, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2620, 0, 3, 2520,
                                                                       1096, 2530, 1192, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2650, 0, 3, 2530,
                                                                       1102, 2540, 1210, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2680, 0, 3, 2540,
                                                                       1108, 2550, 1228, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2710, 0, 3, 2560,
                                                                       1156, 2590, 445, 463,
                                                                       1318, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2770, 0, 3, 2590,
                                                                       1174, 2620, 463, 481,
                                                                       1354, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2830, 0, 3, 2620,
                                                                       1192, 2650, 481, 499,
                                                                       1390, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 2890, 0, 3, 2650,
                                                                       1210, 2680, 499, 517,
                                                                       1426, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 2950, 0, 3, 2710,
                                                                       1318, 2770, 553, 583,
                                                                       1582, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3050, 0, 3, 2770,
                                                                       1354, 2830, 583, 613,
                                                                       1642, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3150, 0, 3, 2830,
                                                                       1390, 2890, 613, 643,
                                                                       1702, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 3250, 0, 3, 2950,
                                                                       1582, 3050, 703, 748,
                                                                       1942, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 3400, 0, 3, 3050,
                                                                       1642, 3150, 748, 793,
                                                                       2032, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 3550, 0, 3, 3250,
                                                                       1942, 3400, 883, 946,
                                                                       2374, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3760, 3, 1072,
                                                                       1078, 2500, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3775, 3, 1078,
                                                                       1084, 2510, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3790, 3, 1084,
                                                                       1090, 2520, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3805, 3, 1090,
                                                                       1096, 2530, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3820, 3, 1096,
                                                                       1102, 2540, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3835, 3, 1102,
                                                                       1108, 2550, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3850, 0, 3, 3760,
                                                                       2500, 3775, 1120, 1138,
                                                                       2560, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3895, 0, 3, 3775,
                                                                       2510, 3790, 1138, 1156,
                                                                       2590, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3940, 0, 3, 3790,
                                                                       2520, 3805, 1156, 1174,
                                                                       2620, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3985, 0, 3, 3805,
                                                                       2530, 3820, 1174, 1192,
                                                                       2650, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4030, 0, 3, 3820,
                                                                       2540, 3835, 1192, 1210,
                                                                       2680, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 4075, 0, 3, 3850,
                                                                       2560, 3895, 1246, 1282,
                                                                       2710, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 4165, 0, 3, 3895,
                                                                       2590, 3940, 1282, 1318,
                                                                       2770, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 4255, 0, 3, 3940,
                                                                       2620, 3985, 1318, 1354,
                                                                       2830, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 4345, 0, 3, 3985,
                                                                       2650, 4030, 1354, 1390,
                                                                       2890, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 4435, 0, 3, 4075,
                                                                       2710, 4165, 1462, 1522,
                                                                       2950, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 4585, 0, 3, 4165,
                                                                       2770, 4255, 1522, 1582,
                                                                       3050, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 4735, 0, 3, 4255,
                                                                       2830, 4345, 1582, 1642,
                                                                       3150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 4885, 0, 3, 4435,
                                                                       2950, 4585, 1762, 1852,
                                                                       3250, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 5110, 0, 3, 4585,
                                                                       3050, 4735, 1852, 1942,
                                                                       3400, ncols, gamma, p,
                                                                       q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 5335, 0, 3, 4885,
                                                                       3250, 5110, 2122, 2248,
                                                                       3550, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 5650, 4885, 225, ncols);

                    simdfunc::contract_primitives(buffer, 6010, 5335, 315, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 5875, 5650, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 6325, 6010, 21, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 6514, 5875, 6325, 9, nmax);

        simdtrf::transform_g_inner(buffer, 6919, 6514, 3, 9, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 6919, 81, nmax);
    }

    for (size_t m = 0; m < 243; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
