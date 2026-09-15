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


#include "SimdThreeCenterElectronRepulsionRecDHD.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_dhd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_dhd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 6178, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 275 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 6178, 3548, 755, dimensions);

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

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 358, 0, 3, 163,
                                                                       178, 253, 274, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 386, 0, 3, 178,
                                                                       193, 274, 295, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 414, 0, 3, 193,
                                                                       208, 295, 316, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 442, 0, 3, 208,
                                                                       223, 316, 337, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 470, 0, 3, 253,
                                                                       274, 358, 386, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 506, 0, 3, 274,
                                                                       295, 386, 414, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 542, 0, 3, 295,
                                                                       316, 414, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 578, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 581, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 584, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 587, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 590, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 593, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 596, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 599, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 602, 3, 10, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 611, 3, 11, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 620, 3, 12, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 629, 3, 13, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 638, 3, 14, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 647, 3, 15, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 656, 3, 16, 42,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 665, 3, 24, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 683, 3, 27, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 701, 3, 30, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 719, 3, 33, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 737, 3, 36, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 755, 3, 39, 87,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 773, 3, 57, 113,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 803, 3, 63, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 833, 3, 69, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 863, 3, 75, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 893, 3, 81, 153,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 923, 3, 113, 193,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 968, 3, 123, 208,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1013, 3, 133, 223,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1058, 3, 143, 238,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1103, 3, 193, 295,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1166, 3, 208, 316,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1229, 3, 223, 337,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1292, 3, 295, 414,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1376, 3, 316, 442,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 1460, 3, 414, 542,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1568, 3, 8, 9,
                                                                       578, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1574, 3, 9, 10,
                                                                       581, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1580, 3, 10, 11,
                                                                       584, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1586, 3, 11, 12,
                                                                       587, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1592, 3, 12, 13,
                                                                       590, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1598, 3, 13, 14,
                                                                       593, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1604, 3, 14, 15,
                                                                       596, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1610, 3, 15, 16,
                                                                       599, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1616, 0, 3, 1568,
                                                                       578, 1574, 602, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1634, 0, 3, 1574,
                                                                       581, 1580, 611, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1652, 0, 3, 1580,
                                                                       584, 1586, 620, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1670, 0, 3, 1586,
                                                                       587, 1592, 629, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 1592,
                                                                       590, 1598, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1706, 0, 3, 1598,
                                                                       593, 1604, 647, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1724, 0, 3, 1604,
                                                                       596, 1610, 656, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1742, 0, 3, 1616,
                                                                       602, 1634, 45, 51, 665,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1778, 0, 3, 1634,
                                                                       611, 1652, 51, 57, 683,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1814, 0, 3, 1652,
                                                                       620, 1670, 57, 63, 701,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1850, 0, 3, 1670,
                                                                       629, 1688, 63, 69, 719,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1886, 0, 3, 1688,
                                                                       638, 1706, 69, 75, 737,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1922, 0, 3, 1706,
                                                                       647, 1724, 75, 81, 755,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1958, 0, 3, 1742,
                                                                       665, 1778, 93, 103, 773,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2018, 0, 3, 1778,
                                                                       683, 1814, 103, 113, 803,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2078, 0, 3, 1814,
                                                                       701, 1850, 113, 123, 833,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2138, 0, 3, 1850,
                                                                       719, 1886, 123, 133, 863,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2198, 0, 3, 1886,
                                                                       737, 1922, 133, 143, 893,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2258, 0, 3, 1958,
                                                                       773, 2018, 163, 178, 923,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2348, 0, 3, 2018,
                                                                       803, 2078, 178, 193, 968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2438, 0, 3, 2078,
                                                                       833, 2138, 193, 208, 1013,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2528, 0, 3, 2138,
                                                                       863, 2198, 208, 223, 1058,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2618, 0, 3, 2258,
                                                                       923, 2348, 253, 274, 1103,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2744, 0, 3, 2348,
                                                                       968, 2438, 274, 295, 1166,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2870, 0, 3, 2438,
                                                                       1013, 2528, 295, 316,
                                                                       1229, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 2996, 0, 3, 2618,
                                                                       1103, 2744, 358, 386,
                                                                       1292, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 3164, 0, 3, 2744,
                                                                       1166, 2870, 386, 414,
                                                                       1376, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 3332, 0, 3, 2996,
                                                                       1292, 3164, 470, 506,
                                                                       1460, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 3548, 2618, 126, ncols);

                    simdfunc::contract_primitives(buffer, 3779, 2996, 168, ncols);

                    simdfunc::contract_primitives(buffer, 4087, 3332, 216, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 3674, 3548, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3947, 3779, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4303, 4087, 36, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 4483, 3674, 3947, 5, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 4798, 3947, 4303, 5, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 5218, 4483, 4798, 5, nmax);

        simdtrf::transform_h_inner(buffer, 5848, 5218, 6, 5, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 5848, 55, nmax);
    }

    for (size_t m = 0; m < 275; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
