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

    const auto nmax = simdfunc::prepare_buffer(buffer, 38813, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 38813, 35357, 1875, dimensions);

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

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 13,
                                                             ncols, fj, mu, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 60, 0, 3, 7, 8,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 66, 0, 3, 8, 9,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 72, 0, 3, 9, 10,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 78, 0, 3, 10, 11,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 84, 0, 3, 11, 12,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 90, 0, 3, 12, 13,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 96, 0, 3, 13, 14,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 102, 0, 3, 14, 15,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 108, 0, 3, 15, 16,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 114, 0, 3, 16, 17,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 120, 0, 3, 17, 18,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 126, 0, 3, 18, 19,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 132, 0, 3, 21, 24,
                                                                       60, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 142, 0, 3, 24, 27,
                                                                       66, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 152, 0, 3, 27, 30,
                                                                       72, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 162, 0, 3, 30, 33,
                                                                       78, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 172, 0, 3, 33, 36,
                                                                       84, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 182, 0, 3, 36, 39,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 39, 42,
                                                                       96, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 202, 0, 3, 42, 45,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 212, 0, 3, 45, 48,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 222, 0, 3, 48, 51,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 232, 0, 3, 51, 54,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 242, 0, 3, 60, 66,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 257, 0, 3, 66, 72,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 272, 0, 3, 72, 78,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 287, 0, 3, 78, 84,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 302, 0, 3, 84, 90,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 317, 0, 3, 90, 96,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 332, 0, 3, 96,
                                                                       102, 192, 202, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 347, 0, 3, 102,
                                                                       108, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 362, 0, 3, 108,
                                                                       114, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 377, 0, 3, 114,
                                                                       120, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 392, 0, 3, 132,
                                                                       142, 242, 257, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 413, 0, 3, 142,
                                                                       152, 257, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 434, 0, 3, 152,
                                                                       162, 272, 287, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 455, 0, 3, 162,
                                                                       172, 287, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 476, 0, 3, 172,
                                                                       182, 302, 317, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 497, 0, 3, 182,
                                                                       192, 317, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 192,
                                                                       202, 332, 347, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 539, 0, 3, 202,
                                                                       212, 347, 362, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 560, 0, 3, 212,
                                                                       222, 362, 377, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 581, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 584, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 587, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 590, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 593, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 596, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 599, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 602, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 605, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 608, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 611, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 614, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 617, 3, 9, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 626, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 635, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 644, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 653, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 662, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 671, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 680, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 689, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 698, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 707, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 716, 3, 27, 72,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 734, 3, 30, 78,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 752, 3, 33, 84,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 770, 3, 36, 90,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 788, 3, 39, 96,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 806, 3, 42, 102,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 824, 3, 45, 108,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 842, 3, 48, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 860, 3, 51, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 878, 3, 54, 126,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 896, 3, 72, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 926, 3, 78, 162,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 956, 3, 84, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 986, 3, 90, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1016, 3, 96, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1046, 3, 102, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1076, 3, 108, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1106, 3, 114, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1136, 3, 120, 232,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1166, 3, 152, 272,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1211, 3, 162, 287,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1256, 3, 172, 302,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1301, 3, 182, 317,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1346, 3, 192, 332,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1391, 3, 202, 347,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1436, 3, 212, 362,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1481, 3, 222, 377,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1526, 3, 272, 434,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1589, 3, 287, 455,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1652, 3, 302, 476,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1715, 3, 317, 497,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1778, 3, 332, 518,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1841, 3, 347, 539,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1904, 3, 362, 560,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1967, 3, 7, 8,
                                                                       581, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1973, 3, 8, 9,
                                                                       584, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1979, 3, 9, 10,
                                                                       587, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1985, 3, 10, 11,
                                                                       590, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1991, 3, 11, 12,
                                                                       593, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1997, 3, 12, 13,
                                                                       596, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2003, 3, 13, 14,
                                                                       599, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2009, 3, 14, 15,
                                                                       602, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2015, 3, 15, 16,
                                                                       605, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2021, 3, 16, 17,
                                                                       608, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2027, 3, 17, 18,
                                                                       611, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2033, 3, 18, 19,
                                                                       614, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2039, 0, 3, 1967,
                                                                       581, 1973, 617, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2057, 0, 3, 1973,
                                                                       584, 1979, 626, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2075, 0, 3, 1979,
                                                                       587, 1985, 635, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2093, 0, 3, 1985,
                                                                       590, 1991, 644, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2111, 0, 3, 1991,
                                                                       593, 1997, 653, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2129, 0, 3, 1997,
                                                                       596, 2003, 662, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2147, 0, 3, 2003,
                                                                       599, 2009, 671, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2165, 0, 3, 2009,
                                                                       602, 2015, 680, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2183, 0, 3, 2015,
                                                                       605, 2021, 689, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2201, 0, 3, 2021,
                                                                       608, 2027, 698, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2219, 0, 3, 2027,
                                                                       611, 2033, 707, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2237, 0, 3, 2039,
                                                                       617, 2057, 60, 66, 716,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2273, 0, 3, 2057,
                                                                       626, 2075, 66, 72, 734,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2309, 0, 3, 2075,
                                                                       635, 2093, 72, 78, 752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2345, 0, 3, 2093,
                                                                       644, 2111, 78, 84, 770,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2381, 0, 3, 2111,
                                                                       653, 2129, 84, 90, 788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2417, 0, 3, 2129,
                                                                       662, 2147, 90, 96, 806,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2453, 0, 3, 2147,
                                                                       671, 2165, 96, 102, 824,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2489, 0, 3, 2165,
                                                                       680, 2183, 102, 108, 842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2525, 0, 3, 2183,
                                                                       689, 2201, 108, 114, 860,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2561, 0, 3, 2201,
                                                                       698, 2219, 114, 120, 878,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2597, 0, 3, 2237,
                                                                       716, 2273, 132, 142, 896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2657, 0, 3, 2273,
                                                                       734, 2309, 142, 152, 926,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2717, 0, 3, 2309,
                                                                       752, 2345, 152, 162, 956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2777, 0, 3, 2345,
                                                                       770, 2381, 162, 172, 986,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2837, 0, 3, 2381,
                                                                       788, 2417, 172, 182, 1016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2897, 0, 3, 2417,
                                                                       806, 2453, 182, 192, 1046,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2957, 0, 3, 2453,
                                                                       824, 2489, 192, 202, 1076,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3017, 0, 3, 2489,
                                                                       842, 2525, 202, 212, 1106,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3077, 0, 3, 2525,
                                                                       860, 2561, 212, 222, 1136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3137, 0, 3, 2597,
                                                                       896, 2657, 242, 257, 1166,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3227, 0, 3, 2657,
                                                                       926, 2717, 257, 272, 1211,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3317, 0, 3, 2717,
                                                                       956, 2777, 272, 287, 1256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3407, 0, 3, 2777,
                                                                       986, 2837, 287, 302, 1301,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3497, 0, 3, 2837,
                                                                       1016, 2897, 302, 317,
                                                                       1346, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3587, 0, 3, 2897,
                                                                       1046, 2957, 317, 332,
                                                                       1391, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3677, 0, 3, 2957,
                                                                       1076, 3017, 332, 347,
                                                                       1436, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3767, 0, 3, 3017,
                                                                       1106, 3077, 347, 362,
                                                                       1481, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3857, 0, 3, 3137,
                                                                       1166, 3227, 392, 413,
                                                                       1526, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3983, 0, 3, 3227,
                                                                       1211, 3317, 413, 434,
                                                                       1589, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4109, 0, 3, 3317,
                                                                       1256, 3407, 434, 455,
                                                                       1652, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4235, 0, 3, 3407,
                                                                       1301, 3497, 455, 476,
                                                                       1715, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4361, 0, 3, 3497,
                                                                       1346, 3587, 476, 497,
                                                                       1778, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4487, 0, 3, 3587,
                                                                       1391, 3677, 497, 518,
                                                                       1841, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4613, 0, 3, 3677,
                                                                       1436, 3767, 518, 539,
                                                                       1904, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4739, 3, 581, 584,
                                                                       1979, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4749, 3, 584, 587,
                                                                       1985, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4759, 3, 587, 590,
                                                                       1991, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4769, 3, 590, 593,
                                                                       1997, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4779, 3, 593, 596,
                                                                       2003, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4789, 3, 596, 599,
                                                                       2009, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4799, 3, 599, 602,
                                                                       2015, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4809, 3, 602, 605,
                                                                       2021, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4819, 3, 605, 608,
                                                                       2027, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4829, 3, 608, 611,
                                                                       2033, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4839, 0, 3, 4739,
                                                                       1979, 4749, 2075, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4869, 0, 3, 4749,
                                                                       1985, 4759, 2093, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4899, 0, 3, 4759,
                                                                       1991, 4769, 2111, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4929, 0, 3, 4769,
                                                                       1997, 4779, 2129, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4959, 0, 3, 4779,
                                                                       2003, 4789, 2147, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4989, 0, 3, 4789,
                                                                       2009, 4799, 2165, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5019, 0, 3, 4799,
                                                                       2015, 4809, 2183, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5049, 0, 3, 4809,
                                                                       2021, 4819, 2201, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5079, 0, 3, 4819,
                                                                       2027, 4829, 2219, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5109, 0, 3, 4839,
                                                                       2075, 4869, 716, 734,
                                                                       2309, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5169, 0, 3, 4869,
                                                                       2093, 4899, 734, 752,
                                                                       2345, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5229, 0, 3, 4899,
                                                                       2111, 4929, 752, 770,
                                                                       2381, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5289, 0, 3, 4929,
                                                                       2129, 4959, 770, 788,
                                                                       2417, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5349, 0, 3, 4959,
                                                                       2147, 4989, 788, 806,
                                                                       2453, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5409, 0, 3, 4989,
                                                                       2165, 5019, 806, 824,
                                                                       2489, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5469, 0, 3, 5019,
                                                                       2183, 5049, 824, 842,
                                                                       2525, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5529, 0, 3, 5049,
                                                                       2201, 5079, 842, 860,
                                                                       2561, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5589, 0, 3, 5109,
                                                                       2309, 5169, 896, 926,
                                                                       2717, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5689, 0, 3, 5169,
                                                                       2345, 5229, 926, 956,
                                                                       2777, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5789, 0, 3, 5229,
                                                                       2381, 5289, 956, 986,
                                                                       2837, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5889, 0, 3, 5289,
                                                                       2417, 5349, 986, 1016,
                                                                       2897, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5989, 0, 3, 5349,
                                                                       2453, 5409, 1016, 1046,
                                                                       2957, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6089, 0, 3, 5409,
                                                                       2489, 5469, 1046, 1076,
                                                                       3017, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6189, 0, 3, 5469,
                                                                       2525, 5529, 1076, 1106,
                                                                       3077, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6289, 0, 3, 5589,
                                                                       2717, 5689, 1166, 1211,
                                                                       3317, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6439, 0, 3, 5689,
                                                                       2777, 5789, 1211, 1256,
                                                                       3407, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6589, 0, 3, 5789,
                                                                       2837, 5889, 1256, 1301,
                                                                       3497, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6739, 0, 3, 5889,
                                                                       2897, 5989, 1301, 1346,
                                                                       3587, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6889, 0, 3, 5989,
                                                                       2957, 6089, 1346, 1391,
                                                                       3677, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7039, 0, 3, 6089,
                                                                       3017, 6189, 1391, 1436,
                                                                       3767, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7189, 0, 3, 6289,
                                                                       3317, 6439, 1526, 1589,
                                                                       4109, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7399, 0, 3, 6439,
                                                                       3407, 6589, 1589, 1652,
                                                                       4235, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7609, 0, 3, 6589,
                                                                       3497, 6739, 1652, 1715,
                                                                       4361, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7819, 0, 3, 6739,
                                                                       3587, 6889, 1715, 1778,
                                                                       4487, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8029, 0, 3, 6889,
                                                                       3677, 7039, 1778, 1841,
                                                                       4613, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8239, 3, 1967,
                                                                       1973, 4739, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8254, 3, 1973,
                                                                       1979, 4749, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8269, 3, 1979,
                                                                       1985, 4759, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8284, 3, 1985,
                                                                       1991, 4769, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8299, 3, 1991,
                                                                       1997, 4779, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8314, 3, 1997,
                                                                       2003, 4789, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8329, 3, 2003,
                                                                       2009, 4799, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8344, 3, 2009,
                                                                       2015, 4809, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8359, 3, 2015,
                                                                       2021, 4819, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8374, 3, 2021,
                                                                       2027, 4829, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8389, 0, 3, 8239,
                                                                       4739, 8254, 2039, 2057,
                                                                       4839, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8434, 0, 3, 8254,
                                                                       4749, 8269, 2057, 2075,
                                                                       4869, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8479, 0, 3, 8269,
                                                                       4759, 8284, 2075, 2093,
                                                                       4899, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8524, 0, 3, 8284,
                                                                       4769, 8299, 2093, 2111,
                                                                       4929, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8569, 0, 3, 8299,
                                                                       4779, 8314, 2111, 2129,
                                                                       4959, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8614, 0, 3, 8314,
                                                                       4789, 8329, 2129, 2147,
                                                                       4989, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8659, 0, 3, 8329,
                                                                       4799, 8344, 2147, 2165,
                                                                       5019, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8704, 0, 3, 8344,
                                                                       4809, 8359, 2165, 2183,
                                                                       5049, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8749, 0, 3, 8359,
                                                                       4819, 8374, 2183, 2201,
                                                                       5079, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8794, 0, 3, 8389,
                                                                       4839, 8434, 2237, 2273,
                                                                       5109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8884, 0, 3, 8434,
                                                                       4869, 8479, 2273, 2309,
                                                                       5169, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8974, 0, 3, 8479,
                                                                       4899, 8524, 2309, 2345,
                                                                       5229, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9064, 0, 3, 8524,
                                                                       4929, 8569, 2345, 2381,
                                                                       5289, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9154, 0, 3, 8569,
                                                                       4959, 8614, 2381, 2417,
                                                                       5349, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9244, 0, 3, 8614,
                                                                       4989, 8659, 2417, 2453,
                                                                       5409, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9334, 0, 3, 8659,
                                                                       5019, 8704, 2453, 2489,
                                                                       5469, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9424, 0, 3, 8704,
                                                                       5049, 8749, 2489, 2525,
                                                                       5529, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9514, 0, 3, 8794,
                                                                       5109, 8884, 2597, 2657,
                                                                       5589, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9664, 0, 3, 8884,
                                                                       5169, 8974, 2657, 2717,
                                                                       5689, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9814, 0, 3, 8974,
                                                                       5229, 9064, 2717, 2777,
                                                                       5789, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9964, 0, 3, 9064,
                                                                       5289, 9154, 2777, 2837,
                                                                       5889, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10114, 0, 3, 9154,
                                                                       5349, 9244, 2837, 2897,
                                                                       5989, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10264, 0, 3, 9244,
                                                                       5409, 9334, 2897, 2957,
                                                                       6089, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10414, 0, 3, 9334,
                                                                       5469, 9424, 2957, 3017,
                                                                       6189, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 10564, 0, 3, 9514,
                                                                       5589, 9664, 3137, 3227,
                                                                       6289, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 10789, 0, 3, 9664,
                                                                       5689, 9814, 3227, 3317,
                                                                       6439, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 11014, 0, 3, 9814,
                                                                       5789, 9964, 3317, 3407,
                                                                       6589, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 11239, 0, 3, 9964,
                                                                       5889, 10114, 3407, 3497,
                                                                       6739, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 11464, 0, 3,
                                                                       10114, 5989, 10264, 3497,
                                                                       3587, 6889, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 11689, 0, 3,
                                                                       10264, 6089, 10414, 3587,
                                                                       3677, 7039, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 11914, 0, 3,
                                                                       10564, 6289, 10789, 3857,
                                                                       3983, 7189, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 12229, 0, 3,
                                                                       10789, 6439, 11014, 3983,
                                                                       4109, 7399, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 12544, 0, 3,
                                                                       11014, 6589, 11239, 4109,
                                                                       4235, 7609, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 12859, 0, 3,
                                                                       11239, 6739, 11464, 4235,
                                                                       4361, 7819, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 13174, 0, 3,
                                                                       11464, 6889, 11689, 4361,
                                                                       4487, 8029, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13489, 3, 4739,
                                                                       4749, 8269, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13510, 3, 4749,
                                                                       4759, 8284, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13531, 3, 4759,
                                                                       4769, 8299, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13552, 3, 4769,
                                                                       4779, 8314, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13573, 3, 4779,
                                                                       4789, 8329, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13594, 3, 4789,
                                                                       4799, 8344, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13615, 3, 4799,
                                                                       4809, 8359, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13636, 3, 4809,
                                                                       4819, 8374, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13657, 0, 3,
                                                                       13489, 8269, 13510, 4839,
                                                                       4869, 8479, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13720, 0, 3,
                                                                       13510, 8284, 13531, 4869,
                                                                       4899, 8524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13783, 0, 3,
                                                                       13531, 8299, 13552, 4899,
                                                                       4929, 8569, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13846, 0, 3,
                                                                       13552, 8314, 13573, 4929,
                                                                       4959, 8614, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13909, 0, 3,
                                                                       13573, 8329, 13594, 4959,
                                                                       4989, 8659, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13972, 0, 3,
                                                                       13594, 8344, 13615, 4989,
                                                                       5019, 8704, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 14035, 0, 3,
                                                                       13615, 8359, 13636, 5019,
                                                                       5049, 8749, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14098, 0, 3,
                                                                       13657, 8479, 13720, 5109,
                                                                       5169, 8974, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14224, 0, 3,
                                                                       13720, 8524, 13783, 5169,
                                                                       5229, 9064, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14350, 0, 3,
                                                                       13783, 8569, 13846, 5229,
                                                                       5289, 9154, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14476, 0, 3,
                                                                       13846, 8614, 13909, 5289,
                                                                       5349, 9244, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14602, 0, 3,
                                                                       13909, 8659, 13972, 5349,
                                                                       5409, 9334, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 14728, 0, 3,
                                                                       13972, 8704, 14035, 5409,
                                                                       5469, 9424, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 14854, 0, 3,
                                                                       14098, 8974, 14224, 5589,
                                                                       5689, 9814, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 15064, 0, 3,
                                                                       14224, 9064, 14350, 5689,
                                                                       5789, 9964, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 15274, 0, 3,
                                                                       14350, 9154, 14476, 5789,
                                                                       5889, 10114, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 15484, 0, 3,
                                                                       14476, 9244, 14602, 5889,
                                                                       5989, 10264, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 15694, 0, 3,
                                                                       14602, 9334, 14728, 5989,
                                                                       6089, 10414, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 15904, 0, 3,
                                                                       14854, 9814, 15064, 6289,
                                                                       6439, 11014, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 16219, 0, 3,
                                                                       15064, 9964, 15274, 6439,
                                                                       6589, 11239, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 16534, 0, 3,
                                                                       15274, 10114, 15484, 6589,
                                                                       6739, 11464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 16849, 0, 3,
                                                                       15484, 10264, 15694, 6739,
                                                                       6889, 11689, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 17164, 0, 3,
                                                                       15904, 11014, 16219, 7189,
                                                                       7399, 12544, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 17605, 0, 3,
                                                                       16219, 11239, 16534, 7399,
                                                                       7609, 12859, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 18046, 0, 3,
                                                                       16534, 11464, 16849, 7609,
                                                                       7819, 13174, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18487, 3, 8239,
                                                                       8254, 13489, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18515, 3, 8254,
                                                                       8269, 13510, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18543, 3, 8269,
                                                                       8284, 13531, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18571, 3, 8284,
                                                                       8299, 13552, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18599, 3, 8299,
                                                                       8314, 13573, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18627, 3, 8314,
                                                                       8329, 13594, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18655, 3, 8329,
                                                                       8344, 13615, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 18683, 3, 8344,
                                                                       8359, 13636, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 18711, 0, 3,
                                                                       18487, 13489, 18515, 8389,
                                                                       8434, 13657, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 18795, 0, 3,
                                                                       18515, 13510, 18543, 8434,
                                                                       8479, 13720, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 18879, 0, 3,
                                                                       18543, 13531, 18571, 8479,
                                                                       8524, 13783, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 18963, 0, 3,
                                                                       18571, 13552, 18599, 8524,
                                                                       8569, 13846, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 19047, 0, 3,
                                                                       18599, 13573, 18627, 8569,
                                                                       8614, 13909, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 19131, 0, 3,
                                                                       18627, 13594, 18655, 8614,
                                                                       8659, 13972, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 19215, 0, 3,
                                                                       18655, 13615, 18683, 8659,
                                                                       8704, 14035, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 19299, 0, 3,
                                                                       18711, 13657, 18795, 8794,
                                                                       8884, 14098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 19467, 0, 3,
                                                                       18795, 13720, 18879, 8884,
                                                                       8974, 14224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 19635, 0, 3,
                                                                       18879, 13783, 18963, 8974,
                                                                       9064, 14350, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 19803, 0, 3,
                                                                       18963, 13846, 19047, 9064,
                                                                       9154, 14476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 19971, 0, 3,
                                                                       19047, 13909, 19131, 9154,
                                                                       9244, 14602, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 20139, 0, 3,
                                                                       19131, 13972, 19215, 9244,
                                                                       9334, 14728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 20307, 0, 3,
                                                                       19299, 14098, 19467, 9514,
                                                                       9664, 14854, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 20587, 0, 3,
                                                                       19467, 14224, 19635, 9664,
                                                                       9814, 15064, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 20867, 0, 3,
                                                                       19635, 14350, 19803, 9814,
                                                                       9964, 15274, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 21147, 0, 3,
                                                                       19803, 14476, 19971, 9964,
                                                                       10114, 15484, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 21427, 0, 3,
                                                                       19971, 14602, 20139,
                                                                       10114, 10264, 15694,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 21707, 0, 3,
                                                                       20307, 14854, 20587,
                                                                       10564, 10789, 15904,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 22127, 0, 3,
                                                                       20587, 15064, 20867,
                                                                       10789, 11014, 16219,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 22547, 0, 3,
                                                                       20867, 15274, 21147,
                                                                       11014, 11239, 16534,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 22967, 0, 3,
                                                                       21147, 15484, 21427,
                                                                       11239, 11464, 16849,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 23387, 0, 3,
                                                                       21707, 15904, 22127,
                                                                       11914, 12229, 17164,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 23975, 0, 3,
                                                                       22127, 16219, 22547,
                                                                       12229, 12544, 17605,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 24563, 0, 3,
                                                                       22547, 16534, 22967,
                                                                       12544, 12859, 18046,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25151, 3, 13489,
                                                                       13510, 18543, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25187, 3, 13510,
                                                                       13531, 18571, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25223, 3, 13531,
                                                                       13552, 18599, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25259, 3, 13552,
                                                                       13573, 18627, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25295, 3, 13573,
                                                                       13594, 18655, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 25331, 3, 13594,
                                                                       13615, 18683, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 25367, 0, 3,
                                                                       25151, 18543, 25187,
                                                                       13657, 13720, 18879,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 25475, 0, 3,
                                                                       25187, 18571, 25223,
                                                                       13720, 13783, 18963,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 25583, 0, 3,
                                                                       25223, 18599, 25259,
                                                                       13783, 13846, 19047,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 25691, 0, 3,
                                                                       25259, 18627, 25295,
                                                                       13846, 13909, 19131,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 25799, 0, 3,
                                                                       25295, 18655, 25331,
                                                                       13909, 13972, 19215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 25907, 0, 3,
                                                                       25367, 18879, 25475,
                                                                       14098, 14224, 19635,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 26123, 0, 3,
                                                                       25475, 18963, 25583,
                                                                       14224, 14350, 19803,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 26339, 0, 3,
                                                                       25583, 19047, 25691,
                                                                       14350, 14476, 19971,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 26555, 0, 3,
                                                                       25691, 19131, 25799,
                                                                       14476, 14602, 20139,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 26771, 0, 3,
                                                                       25907, 19635, 26123,
                                                                       14854, 15064, 20867,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 27131, 0, 3,
                                                                       26123, 19803, 26339,
                                                                       15064, 15274, 21147,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 27491, 0, 3,
                                                                       26339, 19971, 26555,
                                                                       15274, 15484, 21427,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 27851, 0, 3,
                                                                       26771, 20867, 27131,
                                                                       15904, 16219, 22547,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 28391, 0, 3,
                                                                       27131, 21147, 27491,
                                                                       16219, 16534, 22967,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 28931, 0, 3,
                                                                       27851, 22547, 28391,
                                                                       17164, 17605, 24563,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29687, 3, 18487,
                                                                       18515, 25151, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29732, 3, 18515,
                                                                       18543, 25187, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29777, 3, 18543,
                                                                       18571, 25223, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29822, 3, 18571,
                                                                       18599, 25259, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29867, 3, 18599,
                                                                       18627, 25295, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 29912, 3, 18627,
                                                                       18655, 25331, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 29957, 0, 3,
                                                                       29687, 25151, 29732,
                                                                       18711, 18795, 25367,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 30092, 0, 3,
                                                                       29732, 25187, 29777,
                                                                       18795, 18879, 25475,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 30227, 0, 3,
                                                                       29777, 25223, 29822,
                                                                       18879, 18963, 25583,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 30362, 0, 3,
                                                                       29822, 25259, 29867,
                                                                       18963, 19047, 25691,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 30497, 0, 3,
                                                                       29867, 25295, 29912,
                                                                       19047, 19131, 25799,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 30632, 0, 3,
                                                                       29957, 25367, 30092,
                                                                       19299, 19467, 25907,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 30902, 0, 3,
                                                                       30092, 25475, 30227,
                                                                       19467, 19635, 26123,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 31172, 0, 3,
                                                                       30227, 25583, 30362,
                                                                       19635, 19803, 26339,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 31442, 0, 3,
                                                                       30362, 25691, 30497,
                                                                       19803, 19971, 26555,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 31712, 0, 3,
                                                                       30632, 25907, 30902,
                                                                       20307, 20587, 26771,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 32162, 0, 3,
                                                                       30902, 26123, 31172,
                                                                       20587, 20867, 27131,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 32612, 0, 3,
                                                                       31172, 26339, 31442,
                                                                       20867, 21147, 27491,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 33062, 0, 3,
                                                                       31712, 26771, 32162,
                                                                       21707, 22127, 27851,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 33737, 0, 3,
                                                                       32162, 27131, 32612,
                                                                       22127, 22547, 28391,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 34412, 0, 3,
                                                                       33062, 27851, 33737,
                                                                       23387, 23975, 28931,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 35357, 33062, 675, ncols);

                    simdfunc::contract_primitives(buffer, 36287, 34412, 945, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 36032, 35357, 15, 1, nmax);

        simdtrf::transform_l_inner(buffer, 37232, 36287, 21, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 37589, 36032, 37232, 17, nmax);

        simdtrf::transform_g_inner(buffer, 38354, 37589, 3, 17, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 38354, 153, nmax);
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
