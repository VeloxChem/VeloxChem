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


#include "SimdThreeCenterElectronRepulsionRsRecPDL.hpp"

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
#include "SimdTransferPD.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformL.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_pdl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_pdl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 24809, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 510 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 24809, 21958, 1814, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 11,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 19, 3, 11,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 74, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 77, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 80, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 83, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 86, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 89, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 92, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 95, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 98, 0, 3, 7, 8,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 104, 0, 3, 8, 9,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 110, 0, 3, 9, 10,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 116, 0, 3, 10, 11,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 122, 0, 3, 11, 12,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 128, 0, 3, 12, 13,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 134, 0, 3, 13, 14,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 140, 0, 3, 14, 15,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 146, 0, 3, 15, 16,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 152, 0, 3, 16, 17,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 158, 0, 3, 20, 21,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 164, 0, 3, 21, 22,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 170, 0, 3, 22, 23,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 176, 0, 3, 23, 24,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 182, 0, 3, 24, 25,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 188, 0, 3, 25, 26,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 194, 0, 3, 26, 27,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 200, 0, 3, 27, 28,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 206, 0, 3, 28, 29,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 212, 0, 3, 29, 30,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 32, 35,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 35, 38,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 38, 41,
                                                                       110, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 248, 0, 3, 41, 44,
                                                                       116, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 44, 47,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 47, 50,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 50, 53,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 53, 56,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 56, 59,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 65, 68,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 68, 71,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 71, 74,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 74, 77,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 77, 80,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 80, 83,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 83, 86,
                                                                       194, 200, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 86, 89,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 89, 92,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 398, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 401, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 404, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 407, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 410, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 413, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 416, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 419, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 422, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 425, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 428, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 431, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 434, 3, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 437, 3, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 440, 3, 26, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 443, 3, 27, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 446, 3, 28, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 449, 3, 29, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 452, 3, 30, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 455, 3, 31, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 458, 3, 9, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 467, 3, 10, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 476, 3, 11, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 485, 3, 12, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 494, 3, 13, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 503, 3, 14, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 512, 3, 15, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 521, 3, 16, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 530, 3, 17, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 539, 3, 22, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 548, 3, 23, 74,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 557, 3, 24, 77,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 566, 3, 25, 80,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 575, 3, 26, 83,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 584, 3, 27, 86,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 593, 3, 28, 89,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 602, 3, 29, 92,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 611, 3, 30, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 620, 3, 38, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 638, 3, 41, 116,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 656, 3, 44, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 674, 3, 47, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 692, 3, 50, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 710, 3, 53, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 728, 3, 56, 146,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 746, 3, 59, 152,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 764, 3, 71, 170,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 782, 3, 74, 176,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 800, 3, 77, 182,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 818, 3, 80, 188,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 836, 3, 83, 194,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 854, 3, 86, 200,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 872, 3, 89, 206,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 890, 3, 92, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 908, 3, 110, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 938, 3, 116, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 968, 3, 122, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 998, 3, 128, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1028, 3, 134, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1058, 3, 140, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1088, 3, 146, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1118, 3, 170, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1148, 3, 176, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1178, 3, 182, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1208, 3, 188, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1238, 3, 194, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1268, 3, 200, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1298, 3, 206, 388,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1328, 3, 7, 8,
                                                                       398, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1334, 3, 8, 9,
                                                                       401, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1340, 3, 9, 10,
                                                                       404, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1346, 3, 10, 11,
                                                                       407, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1352, 3, 11, 12,
                                                                       410, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1358, 3, 12, 13,
                                                                       413, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1364, 3, 13, 14,
                                                                       416, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1370, 3, 14, 15,
                                                                       419, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1376, 3, 15, 16,
                                                                       422, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1382, 3, 16, 17,
                                                                       425, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1388, 3, 20, 21,
                                                                       428, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1394, 3, 21, 22,
                                                                       431, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1400, 3, 22, 23,
                                                                       434, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1406, 3, 23, 24,
                                                                       437, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1412, 3, 24, 25,
                                                                       440, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1418, 3, 25, 26,
                                                                       443, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1424, 3, 26, 27,
                                                                       446, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1430, 3, 27, 28,
                                                                       449, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1436, 3, 28, 29,
                                                                       452, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1442, 3, 29, 30,
                                                                       455, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1448, 0, 3, 1328,
                                                                       398, 1334, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1466, 0, 3, 1334,
                                                                       401, 1340, 467, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1484, 0, 3, 1340,
                                                                       404, 1346, 476, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1502, 0, 3, 1346,
                                                                       407, 1352, 485, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 1352,
                                                                       410, 1358, 494, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1538, 0, 3, 1358,
                                                                       413, 1364, 503, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1556, 0, 3, 1364,
                                                                       416, 1370, 512, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1574, 0, 3, 1370,
                                                                       419, 1376, 521, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1592, 0, 3, 1376,
                                                                       422, 1382, 530, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1610, 0, 3, 1388,
                                                                       428, 1394, 539, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1628, 0, 3, 1394,
                                                                       431, 1400, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1646, 0, 3, 1400,
                                                                       434, 1406, 557, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1664, 0, 3, 1406,
                                                                       437, 1412, 566, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1682, 0, 3, 1412,
                                                                       440, 1418, 575, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1700, 0, 3, 1418,
                                                                       443, 1424, 584, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1718, 0, 3, 1424,
                                                                       446, 1430, 593, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1736, 0, 3, 1430,
                                                                       449, 1436, 602, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1754, 0, 3, 1436,
                                                                       452, 1442, 611, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 1448,
                                                                       458, 1466, 98, 104, 620,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1808, 0, 3, 1466,
                                                                       467, 1484, 104, 110, 638,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1844, 0, 3, 1484,
                                                                       476, 1502, 110, 116, 656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1880, 0, 3, 1502,
                                                                       485, 1520, 116, 122, 674,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1916, 0, 3, 1520,
                                                                       494, 1538, 122, 128, 692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1952, 0, 3, 1538,
                                                                       503, 1556, 128, 134, 710,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1556,
                                                                       512, 1574, 134, 140, 728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 1574,
                                                                       521, 1592, 140, 146, 746,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2060, 0, 3, 1610,
                                                                       539, 1628, 158, 164, 764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2096, 0, 3, 1628,
                                                                       548, 1646, 164, 170, 782,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2132, 0, 3, 1646,
                                                                       557, 1664, 170, 176, 800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2168, 0, 3, 1664,
                                                                       566, 1682, 176, 182, 818,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2204, 0, 3, 1682,
                                                                       575, 1700, 182, 188, 836,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2240, 0, 3, 1700,
                                                                       584, 1718, 188, 194, 854,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 1718,
                                                                       593, 1736, 194, 200, 872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2312, 0, 3, 1736,
                                                                       602, 1754, 200, 206, 890,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2348, 0, 3, 1772,
                                                                       620, 1808, 218, 228, 908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2408, 0, 3, 1808,
                                                                       638, 1844, 228, 238, 938,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2468, 0, 3, 1844,
                                                                       656, 1880, 238, 248, 968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2528, 0, 3, 1880,
                                                                       674, 1916, 248, 258, 998,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2588, 0, 3, 1916,
                                                                       692, 1952, 258, 268, 1028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1952,
                                                                       710, 1988, 268, 278, 1058,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2708, 0, 3, 1988,
                                                                       728, 2024, 278, 288, 1088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2768, 0, 3, 2060,
                                                                       764, 2096, 308, 318, 1118,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2828, 0, 3, 2096,
                                                                       782, 2132, 318, 328, 1148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2888, 0, 3, 2132,
                                                                       800, 2168, 328, 338, 1178,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2948, 0, 3, 2168,
                                                                       818, 2204, 338, 348, 1208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3008, 0, 3, 2204,
                                                                       836, 2240, 348, 358, 1238,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3068, 0, 3, 2240,
                                                                       854, 2276, 358, 368, 1268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3128, 0, 3, 2276,
                                                                       872, 2312, 368, 378, 1298,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3188, 3, 398, 401,
                                                                       1340, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3198, 3, 401, 404,
                                                                       1346, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3208, 3, 404, 407,
                                                                       1352, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3218, 3, 407, 410,
                                                                       1358, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3228, 3, 410, 413,
                                                                       1364, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3238, 3, 413, 416,
                                                                       1370, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3248, 3, 416, 419,
                                                                       1376, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3258, 3, 419, 422,
                                                                       1382, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3268, 3, 428, 431,
                                                                       1400, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3278, 3, 431, 434,
                                                                       1406, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3288, 3, 434, 437,
                                                                       1412, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3298, 3, 437, 440,
                                                                       1418, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3308, 3, 440, 443,
                                                                       1424, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3318, 3, 443, 446,
                                                                       1430, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3328, 3, 446, 449,
                                                                       1436, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3338, 3, 449, 452,
                                                                       1442, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3348, 0, 3, 3188,
                                                                       1340, 3198, 1484, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3378, 0, 3, 3198,
                                                                       1346, 3208, 1502, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3408, 0, 3, 3208,
                                                                       1352, 3218, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3438, 0, 3, 3218,
                                                                       1358, 3228, 1538, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3468, 0, 3, 3228,
                                                                       1364, 3238, 1556, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3498, 0, 3, 3238,
                                                                       1370, 3248, 1574, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3528, 0, 3, 3248,
                                                                       1376, 3258, 1592, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3558, 0, 3, 3268,
                                                                       1400, 3278, 1646, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3588, 0, 3, 3278,
                                                                       1406, 3288, 1664, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3618, 0, 3, 3288,
                                                                       1412, 3298, 1682, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3648, 0, 3, 3298,
                                                                       1418, 3308, 1700, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3678, 0, 3, 3308,
                                                                       1424, 3318, 1718, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3708, 0, 3, 3318,
                                                                       1430, 3328, 1736, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3738, 0, 3, 3328,
                                                                       1436, 3338, 1754, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3768, 0, 3, 3348,
                                                                       1484, 3378, 620, 638,
                                                                       1844, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3828, 0, 3, 3378,
                                                                       1502, 3408, 638, 656,
                                                                       1880, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3888, 0, 3, 3408,
                                                                       1520, 3438, 656, 674,
                                                                       1916, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3948, 0, 3, 3438,
                                                                       1538, 3468, 674, 692,
                                                                       1952, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4008, 0, 3, 3468,
                                                                       1556, 3498, 692, 710,
                                                                       1988, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4068, 0, 3, 3498,
                                                                       1574, 3528, 710, 728,
                                                                       2024, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4128, 0, 3, 3558,
                                                                       1646, 3588, 764, 782,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4188, 0, 3, 3588,
                                                                       1664, 3618, 782, 800,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4248, 0, 3, 3618,
                                                                       1682, 3648, 800, 818,
                                                                       2204, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4308, 0, 3, 3648,
                                                                       1700, 3678, 818, 836,
                                                                       2240, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4368, 0, 3, 3678,
                                                                       1718, 3708, 836, 854,
                                                                       2276, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4428, 0, 3, 3708,
                                                                       1736, 3738, 854, 872,
                                                                       2312, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4488, 0, 3, 3768,
                                                                       1844, 3828, 908, 938,
                                                                       2468, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4588, 0, 3, 3828,
                                                                       1880, 3888, 938, 968,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4688, 0, 3, 3888,
                                                                       1916, 3948, 968, 998,
                                                                       2588, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4788, 0, 3, 3948,
                                                                       1952, 4008, 998, 1028,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4888, 0, 3, 4008,
                                                                       1988, 4068, 1028, 1058,
                                                                       2708, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4988, 0, 3, 4128,
                                                                       2132, 4188, 1118, 1148,
                                                                       2888, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5088, 0, 3, 4188,
                                                                       2168, 4248, 1148, 1178,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5188, 0, 3, 4248,
                                                                       2204, 4308, 1178, 1208,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 4308,
                                                                       2240, 4368, 1208, 1238,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5388, 0, 3, 4368,
                                                                       2276, 4428, 1238, 1268,
                                                                       3128, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5488, 3, 1328,
                                                                       1334, 3188, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5503, 3, 1334,
                                                                       1340, 3198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5518, 3, 1340,
                                                                       1346, 3208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5533, 3, 1346,
                                                                       1352, 3218, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5548, 3, 1352,
                                                                       1358, 3228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5563, 3, 1358,
                                                                       1364, 3238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5578, 3, 1364,
                                                                       1370, 3248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5593, 3, 1370,
                                                                       1376, 3258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5608, 3, 1388,
                                                                       1394, 3268, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5623, 3, 1394,
                                                                       1400, 3278, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5638, 3, 1400,
                                                                       1406, 3288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5653, 3, 1406,
                                                                       1412, 3298, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5668, 3, 1412,
                                                                       1418, 3308, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5683, 3, 1418,
                                                                       1424, 3318, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5698, 3, 1424,
                                                                       1430, 3328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5713, 3, 1430,
                                                                       1436, 3338, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5728, 0, 3, 5488,
                                                                       3188, 5503, 1448, 1466,
                                                                       3348, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5773, 0, 3, 5503,
                                                                       3198, 5518, 1466, 1484,
                                                                       3378, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5818, 0, 3, 5518,
                                                                       3208, 5533, 1484, 1502,
                                                                       3408, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5863, 0, 3, 5533,
                                                                       3218, 5548, 1502, 1520,
                                                                       3438, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5908, 0, 3, 5548,
                                                                       3228, 5563, 1520, 1538,
                                                                       3468, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5953, 0, 3, 5563,
                                                                       3238, 5578, 1538, 1556,
                                                                       3498, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5998, 0, 3, 5578,
                                                                       3248, 5593, 1556, 1574,
                                                                       3528, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6043, 0, 3, 5608,
                                                                       3268, 5623, 1610, 1628,
                                                                       3558, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6088, 0, 3, 5623,
                                                                       3278, 5638, 1628, 1646,
                                                                       3588, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6133, 0, 3, 5638,
                                                                       3288, 5653, 1646, 1664,
                                                                       3618, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6178, 0, 3, 5653,
                                                                       3298, 5668, 1664, 1682,
                                                                       3648, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6223, 0, 3, 5668,
                                                                       3308, 5683, 1682, 1700,
                                                                       3678, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6268, 0, 3, 5683,
                                                                       3318, 5698, 1700, 1718,
                                                                       3708, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6313, 0, 3, 5698,
                                                                       3328, 5713, 1718, 1736,
                                                                       3738, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6358, 0, 3, 5728,
                                                                       3348, 5773, 1772, 1808,
                                                                       3768, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6448, 0, 3, 5773,
                                                                       3378, 5818, 1808, 1844,
                                                                       3828, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6538, 0, 3, 5818,
                                                                       3408, 5863, 1844, 1880,
                                                                       3888, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6628, 0, 3, 5863,
                                                                       3438, 5908, 1880, 1916,
                                                                       3948, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6718, 0, 3, 5908,
                                                                       3468, 5953, 1916, 1952,
                                                                       4008, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6808, 0, 3, 5953,
                                                                       3498, 5998, 1952, 1988,
                                                                       4068, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6898, 0, 3, 6043,
                                                                       3558, 6088, 2060, 2096,
                                                                       4128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 6988, 0, 3, 6088,
                                                                       3588, 6133, 2096, 2132,
                                                                       4188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7078, 0, 3, 6133,
                                                                       3618, 6178, 2132, 2168,
                                                                       4248, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7168, 0, 3, 6178,
                                                                       3648, 6223, 2168, 2204,
                                                                       4308, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7258, 0, 3, 6223,
                                                                       3678, 6268, 2204, 2240,
                                                                       4368, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7348, 0, 3, 6268,
                                                                       3708, 6313, 2240, 2276,
                                                                       4428, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 7438, 0, 3, 6358,
                                                                       3768, 6448, 2348, 2408,
                                                                       4488, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 7588, 0, 3, 6448,
                                                                       3828, 6538, 2408, 2468,
                                                                       4588, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 7738, 0, 3, 6538,
                                                                       3888, 6628, 2468, 2528,
                                                                       4688, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 7888, 0, 3, 6628,
                                                                       3948, 6718, 2528, 2588,
                                                                       4788, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8038, 0, 3, 6718,
                                                                       4008, 6808, 2588, 2648,
                                                                       4888, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8188, 0, 3, 6898,
                                                                       4128, 6988, 2768, 2828,
                                                                       4988, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8338, 0, 3, 6988,
                                                                       4188, 7078, 2828, 2888,
                                                                       5088, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8488, 0, 3, 7078,
                                                                       4248, 7168, 2888, 2948,
                                                                       5188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8638, 0, 3, 7168,
                                                                       4308, 7258, 2948, 3008,
                                                                       5288, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8788, 0, 3, 7258,
                                                                       4368, 7348, 3008, 3068,
                                                                       5388, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8938, 3, 3188,
                                                                       3198, 5518, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8959, 3, 3198,
                                                                       3208, 5533, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8980, 3, 3208,
                                                                       3218, 5548, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9001, 3, 3218,
                                                                       3228, 5563, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9022, 3, 3228,
                                                                       3238, 5578, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9043, 3, 3238,
                                                                       3248, 5593, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9064, 3, 3268,
                                                                       3278, 5638, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9085, 3, 3278,
                                                                       3288, 5653, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9106, 3, 3288,
                                                                       3298, 5668, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9127, 3, 3298,
                                                                       3308, 5683, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9148, 3, 3308,
                                                                       3318, 5698, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9169, 3, 3318,
                                                                       3328, 5713, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9190, 0, 3, 8938,
                                                                       5518, 8959, 3348, 3378,
                                                                       5818, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9253, 0, 3, 8959,
                                                                       5533, 8980, 3378, 3408,
                                                                       5863, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9316, 0, 3, 8980,
                                                                       5548, 9001, 3408, 3438,
                                                                       5908, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9379, 0, 3, 9001,
                                                                       5563, 9022, 3438, 3468,
                                                                       5953, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9442, 0, 3, 9022,
                                                                       5578, 9043, 3468, 3498,
                                                                       5998, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9505, 0, 3, 9064,
                                                                       5638, 9085, 3558, 3588,
                                                                       6133, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9568, 0, 3, 9085,
                                                                       5653, 9106, 3588, 3618,
                                                                       6178, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9631, 0, 3, 9106,
                                                                       5668, 9127, 3618, 3648,
                                                                       6223, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9694, 0, 3, 9127,
                                                                       5683, 9148, 3648, 3678,
                                                                       6268, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9757, 0, 3, 9148,
                                                                       5698, 9169, 3678, 3708,
                                                                       6313, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9820, 0, 3, 9190,
                                                                       5818, 9253, 3768, 3828,
                                                                       6538, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9946, 0, 3, 9253,
                                                                       5863, 9316, 3828, 3888,
                                                                       6628, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 10072, 0, 3, 9316,
                                                                       5908, 9379, 3888, 3948,
                                                                       6718, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 10198, 0, 3, 9379,
                                                                       5953, 9442, 3948, 4008,
                                                                       6808, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 10324, 0, 3, 9505,
                                                                       6133, 9568, 4128, 4188,
                                                                       7078, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 10450, 0, 3, 9568,
                                                                       6178, 9631, 4188, 4248,
                                                                       7168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 10576, 0, 3, 9631,
                                                                       6223, 9694, 4248, 4308,
                                                                       7258, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 10702, 0, 3, 9694,
                                                                       6268, 9757, 4308, 4368,
                                                                       7348, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 10828, 0, 3, 9820,
                                                                       6538, 9946, 4488, 4588,
                                                                       7738, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 11038, 0, 3, 9946,
                                                                       6628, 10072, 4588, 4688,
                                                                       7888, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 11248, 0, 3,
                                                                       10072, 6718, 10198, 4688,
                                                                       4788, 8038, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 11458, 0, 3,
                                                                       10324, 7078, 10450, 4988,
                                                                       5088, 8488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 11668, 0, 3,
                                                                       10450, 7168, 10576, 5088,
                                                                       5188, 8638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 11878, 0, 3,
                                                                       10576, 7258, 10702, 5188,
                                                                       5288, 8788, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12088, 3, 5488,
                                                                       5503, 8938, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12116, 3, 5503,
                                                                       5518, 8959, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12144, 3, 5518,
                                                                       5533, 8980, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12172, 3, 5533,
                                                                       5548, 9001, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12200, 3, 5548,
                                                                       5563, 9022, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12228, 3, 5563,
                                                                       5578, 9043, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12256, 3, 5608,
                                                                       5623, 9064, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12284, 3, 5623,
                                                                       5638, 9085, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12312, 3, 5638,
                                                                       5653, 9106, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12340, 3, 5653,
                                                                       5668, 9127, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12368, 3, 5668,
                                                                       5683, 9148, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 12396, 3, 5683,
                                                                       5698, 9169, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 12424, 0, 3,
                                                                       12088, 8938, 12116, 5728,
                                                                       5773, 9190, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 12508, 0, 3,
                                                                       12116, 8959, 12144, 5773,
                                                                       5818, 9253, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 12592, 0, 3,
                                                                       12144, 8980, 12172, 5818,
                                                                       5863, 9316, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 12676, 0, 3,
                                                                       12172, 9001, 12200, 5863,
                                                                       5908, 9379, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 12760, 0, 3,
                                                                       12200, 9022, 12228, 5908,
                                                                       5953, 9442, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 12844, 0, 3,
                                                                       12256, 9064, 12284, 6043,
                                                                       6088, 9505, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 12928, 0, 3,
                                                                       12284, 9085, 12312, 6088,
                                                                       6133, 9568, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 13012, 0, 3,
                                                                       12312, 9106, 12340, 6133,
                                                                       6178, 9631, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 13096, 0, 3,
                                                                       12340, 9127, 12368, 6178,
                                                                       6223, 9694, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 13180, 0, 3,
                                                                       12368, 9148, 12396, 6223,
                                                                       6268, 9757, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 13264, 0, 3,
                                                                       12424, 9190, 12508, 6358,
                                                                       6448, 9820, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 13432, 0, 3,
                                                                       12508, 9253, 12592, 6448,
                                                                       6538, 9946, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 13600, 0, 3,
                                                                       12592, 9316, 12676, 6538,
                                                                       6628, 10072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 13768, 0, 3,
                                                                       12676, 9379, 12760, 6628,
                                                                       6718, 10198, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 13936, 0, 3,
                                                                       12844, 9505, 12928, 6898,
                                                                       6988, 10324, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 14104, 0, 3,
                                                                       12928, 9568, 13012, 6988,
                                                                       7078, 10450, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 14272, 0, 3,
                                                                       13012, 9631, 13096, 7078,
                                                                       7168, 10576, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 14440, 0, 3,
                                                                       13096, 9694, 13180, 7168,
                                                                       7258, 10702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 14608, 0, 3,
                                                                       13264, 9820, 13432, 7438,
                                                                       7588, 10828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 14888, 0, 3,
                                                                       13432, 9946, 13600, 7588,
                                                                       7738, 11038, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 15168, 0, 3,
                                                                       13600, 10072, 13768, 7738,
                                                                       7888, 11248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 15448, 0, 3,
                                                                       13936, 10324, 14104, 8188,
                                                                       8338, 11458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 15728, 0, 3,
                                                                       14104, 10450, 14272, 8338,
                                                                       8488, 11668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 16008, 0, 3,
                                                                       14272, 10576, 14440, 8488,
                                                                       8638, 11878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 16288, 3, 8938,
                                                                       8959, 12144, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 16324, 3, 8959,
                                                                       8980, 12172, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 16360, 3, 8980,
                                                                       9001, 12200, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 16396, 3, 9001,
                                                                       9022, 12228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 16432, 3, 9064,
                                                                       9085, 12312, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 16468, 3, 9085,
                                                                       9106, 12340, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 16504, 3, 9106,
                                                                       9127, 12368, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 16540, 3, 9127,
                                                                       9148, 12396, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 16576, 0, 3,
                                                                       16288, 12144, 16324, 9190,
                                                                       9253, 12592, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 16684, 0, 3,
                                                                       16324, 12172, 16360, 9253,
                                                                       9316, 12676, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 16792, 0, 3,
                                                                       16360, 12200, 16396, 9316,
                                                                       9379, 12760, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 16900, 0, 3,
                                                                       16432, 12312, 16468, 9505,
                                                                       9568, 13012, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 17008, 0, 3,
                                                                       16468, 12340, 16504, 9568,
                                                                       9631, 13096, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 17116, 0, 3,
                                                                       16504, 12368, 16540, 9631,
                                                                       9694, 13180, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 17224, 0, 3,
                                                                       16576, 12592, 16684, 9820,
                                                                       9946, 13600, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 17440, 0, 3,
                                                                       16684, 12676, 16792, 9946,
                                                                       10072, 13768, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 17656, 0, 3,
                                                                       16900, 13012, 17008,
                                                                       10324, 10450, 14272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 17872, 0, 3,
                                                                       17008, 13096, 17116,
                                                                       10450, 10576, 14440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 18088, 0, 3,
                                                                       17224, 13600, 17440,
                                                                       10828, 11038, 15168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 18448, 0, 3,
                                                                       17656, 14272, 17872,
                                                                       11458, 11668, 16008,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 18808, 3, 12088,
                                                                       12116, 16288, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 18853, 3, 12116,
                                                                       12144, 16324, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 18898, 3, 12144,
                                                                       12172, 16360, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 18943, 3, 12172,
                                                                       12200, 16396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 18988, 3, 12256,
                                                                       12284, 16432, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 19033, 3, 12284,
                                                                       12312, 16468, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 19078, 3, 12312,
                                                                       12340, 16504, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 19123, 3, 12340,
                                                                       12368, 16540, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 19168, 0, 3,
                                                                       18808, 16288, 18853,
                                                                       12424, 12508, 16576,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 19303, 0, 3,
                                                                       18853, 16324, 18898,
                                                                       12508, 12592, 16684,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 19438, 0, 3,
                                                                       18898, 16360, 18943,
                                                                       12592, 12676, 16792,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 19573, 0, 3,
                                                                       18988, 16432, 19033,
                                                                       12844, 12928, 16900,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 19708, 0, 3,
                                                                       19033, 16468, 19078,
                                                                       12928, 13012, 17008,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 19843, 0, 3,
                                                                       19078, 16504, 19123,
                                                                       13012, 13096, 17116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 19978, 0, 3,
                                                                       19168, 16576, 19303,
                                                                       13264, 13432, 17224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 20248, 0, 3,
                                                                       19303, 16684, 19438,
                                                                       13432, 13600, 17440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 20518, 0, 3,
                                                                       19573, 16900, 19708,
                                                                       13936, 14104, 17656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 20788, 0, 3,
                                                                       19708, 17008, 19843,
                                                                       14104, 14272, 17872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 21058, 0, 3,
                                                                       19978, 17224, 20248,
                                                                       14608, 14888, 18088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 21508, 0, 3,
                                                                       20518, 17656, 20788,
                                                                       15448, 15728, 18448,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 21958, 19978, 270, ncols);

                    simdfunc::contract_primitives(buffer, 22330, 20518, 270, ncols);

                    simdfunc::contract_primitives(buffer, 22702, 21058, 450, ncols);

                    simdfunc::contract_primitives(buffer, 23322, 21508, 450, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 22228, 21958, 6, 1, nmax);

        simdtrf::transform_l_inner(buffer, 22600, 22330, 6, 1, nmax);

        simdtrf::transform_l_inner(buffer, 23152, 22702, 10, 1, nmax);

        simdtrf::transform_l_inner(buffer, 23772, 23322, 10, 1, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 23942, 22228, 23152, 17, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 24248, 22600, 23772, 17, nmax);

        simdtrf::transform_d_inner(buffer, 24554, 24248, 3, 17, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 24554, 85, nmax);

        simdtrf::transform_d_inner(buffer, 24554, 23942, 3, 17, nmax);

        simdtrf::transform_p_outer(values + 255 * nvalues + n * npairs, nvalues, buffer, 24554,
                                   85, nmax);
    }

    for (size_t m = 0; m < 510; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
