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


#include "SimdThreeCenterElectronRepulsionRecPGK.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
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
#include "SimdTransferPG.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_pgk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_pgk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 26986, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 405 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 26986, 24070, 1521, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12}, ncols, fj, 6,
                                                        fq);

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

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 53, 59,
                                                                       113, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 59, 65,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 233, 0, 3, 65, 71,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 248, 0, 3, 71, 77,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 263, 0, 3, 77, 83,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 83, 89,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 293, 0, 3, 89, 95,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 95,
                                                                       101, 183, 193, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 323, 0, 3, 113,
                                                                       123, 203, 218, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 344, 0, 3, 123,
                                                                       133, 218, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 365, 0, 3, 133,
                                                                       143, 233, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 386, 0, 3, 143,
                                                                       153, 248, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 407, 0, 3, 153,
                                                                       163, 263, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 163,
                                                                       173, 278, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 449, 0, 3, 173,
                                                                       183, 293, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 470, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 473, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 476, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 479, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 482, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 485, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 488, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 491, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 494, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 497, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 500, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 503, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 506, 3, 10, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 515, 3, 11, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 524, 3, 12, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 533, 3, 13, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 542, 3, 14, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 551, 3, 15, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 560, 3, 16, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 569, 3, 17, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 578, 3, 18, 50,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 587, 3, 20, 53,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 605, 3, 23, 59,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 623, 3, 26, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 641, 3, 29, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 659, 3, 32, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 677, 3, 35, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 695, 3, 38, 89,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 713, 3, 41, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 731, 3, 44, 101,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 749, 3, 47, 107,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 767, 3, 53, 113,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 797, 3, 59, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 827, 3, 65, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 857, 3, 71, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 887, 3, 77, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 917, 3, 83, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 947, 3, 89, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 977, 3, 95, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1007, 3, 101, 193,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1037, 3, 113, 203,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1082, 3, 123, 218,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1127, 3, 133, 233,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1172, 3, 143, 248,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1217, 3, 153, 263,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1262, 3, 163, 278,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1307, 3, 173, 293,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1352, 3, 183, 308,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1397, 3, 203, 323,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1460, 3, 218, 344,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1523, 3, 233, 365,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1586, 3, 248, 386,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1649, 3, 263, 407,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1712, 3, 278, 428,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1775, 3, 293, 449,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1838, 3, 8, 9,
                                                                       476, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1844, 3, 9, 10,
                                                                       479, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1850, 3, 10, 11,
                                                                       482, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1856, 3, 11, 12,
                                                                       485, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1862, 3, 12, 13,
                                                                       488, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1868, 3, 13, 14,
                                                                       491, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1874, 3, 14, 15,
                                                                       494, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1880, 3, 15, 16,
                                                                       497, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1886, 3, 16, 17,
                                                                       500, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1892, 3, 17, 18,
                                                                       503, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1898, 0, 3, 1838,
                                                                       476, 1844, 506, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1916, 0, 3, 1844,
                                                                       479, 1850, 515, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1934, 0, 3, 1850,
                                                                       482, 1856, 524, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1952, 0, 3, 1856,
                                                                       485, 1862, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1970, 0, 3, 1862,
                                                                       488, 1868, 542, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1868,
                                                                       491, 1874, 551, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2006, 0, 3, 1874,
                                                                       494, 1880, 560, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 1880,
                                                                       497, 1886, 569, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2042, 0, 3, 1886,
                                                                       500, 1892, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2060, 0, 3, 1898,
                                                                       506, 1916, 53, 59, 623,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2096, 0, 3, 1916,
                                                                       515, 1934, 59, 65, 641,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2132, 0, 3, 1934,
                                                                       524, 1952, 65, 71, 659,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2168, 0, 3, 1952,
                                                                       533, 1970, 71, 77, 677,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2204, 0, 3, 1970,
                                                                       542, 1988, 77, 83, 695,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2240, 0, 3, 1988,
                                                                       551, 2006, 83, 89, 713,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 2006,
                                                                       560, 2024, 89, 95, 731,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2312, 0, 3, 2024,
                                                                       569, 2042, 95, 101, 749,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2348, 0, 3, 2060,
                                                                       623, 2096, 113, 123, 827,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2408, 0, 3, 2096,
                                                                       641, 2132, 123, 133, 857,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2468, 0, 3, 2132,
                                                                       659, 2168, 133, 143, 887,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2528, 0, 3, 2168,
                                                                       677, 2204, 143, 153, 917,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2588, 0, 3, 2204,
                                                                       695, 2240, 153, 163, 947,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 2240,
                                                                       713, 2276, 163, 173, 977,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2708, 0, 3, 2276,
                                                                       731, 2312, 173, 183, 1007,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2768, 0, 3, 2348,
                                                                       827, 2408, 203, 218, 1127,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2858, 0, 3, 2408,
                                                                       857, 2468, 218, 233, 1172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2948, 0, 3, 2468,
                                                                       887, 2528, 233, 248, 1217,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3038, 0, 3, 2528,
                                                                       917, 2588, 248, 263, 1262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3128, 0, 3, 2588,
                                                                       947, 2648, 263, 278, 1307,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3218, 0, 3, 2648,
                                                                       977, 2708, 278, 293, 1352,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 2768,
                                                                       1127, 2858, 323, 344,
                                                                       1523, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3434, 0, 3, 2858,
                                                                       1172, 2948, 344, 365,
                                                                       1586, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3560, 0, 3, 2948,
                                                                       1217, 3038, 365, 386,
                                                                       1649, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3686, 0, 3, 3038,
                                                                       1262, 3128, 386, 407,
                                                                       1712, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3812, 0, 3, 3128,
                                                                       1307, 3218, 407, 428,
                                                                       1775, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3938, 3, 470, 473,
                                                                       1838, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3948, 3, 473, 476,
                                                                       1844, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3958, 3, 476, 479,
                                                                       1850, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3968, 3, 479, 482,
                                                                       1856, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3978, 3, 482, 485,
                                                                       1862, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3988, 3, 485, 488,
                                                                       1868, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3998, 3, 488, 491,
                                                                       1874, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4008, 3, 491, 494,
                                                                       1880, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4018, 3, 494, 497,
                                                                       1886, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4028, 3, 497, 500,
                                                                       1892, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4038, 0, 3, 3938,
                                                                       1838, 3948, 1898, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4068, 0, 3, 3948,
                                                                       1844, 3958, 1916, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4098, 0, 3, 3958,
                                                                       1850, 3968, 1934, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4128, 0, 3, 3968,
                                                                       1856, 3978, 1952, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4158, 0, 3, 3978,
                                                                       1862, 3988, 1970, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4188, 0, 3, 3988,
                                                                       1868, 3998, 1988, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4218, 0, 3, 3998,
                                                                       1874, 4008, 2006, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4248, 0, 3, 4008,
                                                                       1880, 4018, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4278, 0, 3, 4018,
                                                                       1886, 4028, 2042, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4308, 0, 3, 4038,
                                                                       1898, 4068, 587, 605,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4368, 0, 3, 4068,
                                                                       1916, 4098, 605, 623,
                                                                       2096, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4428, 0, 3, 4098,
                                                                       1934, 4128, 623, 641,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4488, 0, 3, 4128,
                                                                       1952, 4158, 641, 659,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4548, 0, 3, 4158,
                                                                       1970, 4188, 659, 677,
                                                                       2204, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4608, 0, 3, 4188,
                                                                       1988, 4218, 677, 695,
                                                                       2240, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4668, 0, 3, 4218,
                                                                       2006, 4248, 695, 713,
                                                                       2276, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4728, 0, 3, 4248,
                                                                       2024, 4278, 713, 731,
                                                                       2312, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4788, 0, 3, 4308,
                                                                       2060, 4368, 767, 797,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4888, 0, 3, 4368,
                                                                       2096, 4428, 797, 827,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4988, 0, 3, 4428,
                                                                       2132, 4488, 827, 857,
                                                                       2468, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5088, 0, 3, 4488,
                                                                       2168, 4548, 857, 887,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5188, 0, 3, 4548,
                                                                       2204, 4608, 887, 917,
                                                                       2588, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 4608,
                                                                       2240, 4668, 917, 947,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5388, 0, 3, 4668,
                                                                       2276, 4728, 947, 977,
                                                                       2708, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5488, 0, 3, 4788,
                                                                       2348, 4888, 1037, 1082,
                                                                       2768, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5638, 0, 3, 4888,
                                                                       2408, 4988, 1082, 1127,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5788, 0, 3, 4988,
                                                                       2468, 5088, 1127, 1172,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5938, 0, 3, 5088,
                                                                       2528, 5188, 1172, 1217,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6088, 0, 3, 5188,
                                                                       2588, 5288, 1217, 1262,
                                                                       3128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 6238, 0, 3, 5288,
                                                                       2648, 5388, 1262, 1307,
                                                                       3218, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 6388, 0, 3, 5488,
                                                                       2768, 5638, 1397, 1460,
                                                                       3308, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 6598, 0, 3, 5638,
                                                                       2858, 5788, 1460, 1523,
                                                                       3434, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 6808, 0, 3, 5788,
                                                                       2948, 5938, 1523, 1586,
                                                                       3560, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7018, 0, 3, 5938,
                                                                       3038, 6088, 1586, 1649,
                                                                       3686, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7228, 0, 3, 6088,
                                                                       3128, 6238, 1649, 1712,
                                                                       3812, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7438, 3, 1838,
                                                                       1844, 3958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7453, 3, 1844,
                                                                       1850, 3968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7468, 3, 1850,
                                                                       1856, 3978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7483, 3, 1856,
                                                                       1862, 3988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7498, 3, 1862,
                                                                       1868, 3998, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7513, 3, 1868,
                                                                       1874, 4008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7528, 3, 1874,
                                                                       1880, 4018, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7543, 3, 1880,
                                                                       1886, 4028, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7558, 0, 3, 7438,
                                                                       3958, 7453, 1898, 1916,
                                                                       4098, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7603, 0, 3, 7453,
                                                                       3968, 7468, 1916, 1934,
                                                                       4128, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7648, 0, 3, 7468,
                                                                       3978, 7483, 1934, 1952,
                                                                       4158, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7693, 0, 3, 7483,
                                                                       3988, 7498, 1952, 1970,
                                                                       4188, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7738, 0, 3, 7498,
                                                                       3998, 7513, 1970, 1988,
                                                                       4218, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7783, 0, 3, 7513,
                                                                       4008, 7528, 1988, 2006,
                                                                       4248, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7828, 0, 3, 7528,
                                                                       4018, 7543, 2006, 2024,
                                                                       4278, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7873, 0, 3, 7558,
                                                                       4098, 7603, 2060, 2096,
                                                                       4428, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7963, 0, 3, 7603,
                                                                       4128, 7648, 2096, 2132,
                                                                       4488, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8053, 0, 3, 7648,
                                                                       4158, 7693, 2132, 2168,
                                                                       4548, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8143, 0, 3, 7693,
                                                                       4188, 7738, 2168, 2204,
                                                                       4608, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8233, 0, 3, 7738,
                                                                       4218, 7783, 2204, 2240,
                                                                       4668, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8323, 0, 3, 7783,
                                                                       4248, 7828, 2240, 2276,
                                                                       4728, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8413, 0, 3, 7873,
                                                                       4428, 7963, 2348, 2408,
                                                                       4988, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8563, 0, 3, 7963,
                                                                       4488, 8053, 2408, 2468,
                                                                       5088, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8713, 0, 3, 8053,
                                                                       4548, 8143, 2468, 2528,
                                                                       5188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 8863, 0, 3, 8143,
                                                                       4608, 8233, 2528, 2588,
                                                                       5288, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 9013, 0, 3, 8233,
                                                                       4668, 8323, 2588, 2648,
                                                                       5388, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 9163, 0, 3, 8413,
                                                                       4988, 8563, 2768, 2858,
                                                                       5788, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 9388, 0, 3, 8563,
                                                                       5088, 8713, 2858, 2948,
                                                                       5938, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 9613, 0, 3, 8713,
                                                                       5188, 8863, 2948, 3038,
                                                                       6088, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 9838, 0, 3, 8863,
                                                                       5288, 9013, 3038, 3128,
                                                                       6238, ncols, gamma, p,
                                                                       q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 10063, 0, 3, 9163,
                                                                       5788, 9388, 3308, 3434,
                                                                       6808, ncols, gamma, p,
                                                                       q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 10378, 0, 3, 9388,
                                                                       5938, 9613, 3434, 3560,
                                                                       7018, ncols, gamma, p,
                                                                       q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 10693, 0, 3, 9613,
                                                                       6088, 9838, 3560, 3686,
                                                                       7228, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11008, 3, 3938,
                                                                       3948, 7438, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11029, 3, 3948,
                                                                       3958, 7453, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11050, 3, 3958,
                                                                       3968, 7468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11071, 3, 3968,
                                                                       3978, 7483, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11092, 3, 3978,
                                                                       3988, 7498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11113, 3, 3988,
                                                                       3998, 7513, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11134, 3, 3998,
                                                                       4008, 7528, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11155, 3, 4008,
                                                                       4018, 7543, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11176, 0, 3,
                                                                       11008, 7438, 11029, 4038,
                                                                       4068, 7558, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11239, 0, 3,
                                                                       11029, 7453, 11050, 4068,
                                                                       4098, 7603, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11302, 0, 3,
                                                                       11050, 7468, 11071, 4098,
                                                                       4128, 7648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11365, 0, 3,
                                                                       11071, 7483, 11092, 4128,
                                                                       4158, 7693, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11428, 0, 3,
                                                                       11092, 7498, 11113, 4158,
                                                                       4188, 7738, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11491, 0, 3,
                                                                       11113, 7513, 11134, 4188,
                                                                       4218, 7783, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 11554, 0, 3,
                                                                       11134, 7528, 11155, 4218,
                                                                       4248, 7828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 11617, 0, 3,
                                                                       11176, 7558, 11239, 4308,
                                                                       4368, 7873, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 11743, 0, 3,
                                                                       11239, 7603, 11302, 4368,
                                                                       4428, 7963, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 11869, 0, 3,
                                                                       11302, 7648, 11365, 4428,
                                                                       4488, 8053, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 11995, 0, 3,
                                                                       11365, 7693, 11428, 4488,
                                                                       4548, 8143, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 12121, 0, 3,
                                                                       11428, 7738, 11491, 4548,
                                                                       4608, 8233, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 12247, 0, 3,
                                                                       11491, 7783, 11554, 4608,
                                                                       4668, 8323, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 12373, 0, 3,
                                                                       11617, 7873, 11743, 4788,
                                                                       4888, 8413, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 12583, 0, 3,
                                                                       11743, 7963, 11869, 4888,
                                                                       4988, 8563, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 12793, 0, 3,
                                                                       11869, 8053, 11995, 4988,
                                                                       5088, 8713, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 13003, 0, 3,
                                                                       11995, 8143, 12121, 5088,
                                                                       5188, 8863, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 13213, 0, 3,
                                                                       12121, 8233, 12247, 5188,
                                                                       5288, 9013, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 13423, 0, 3,
                                                                       12373, 8413, 12583, 5488,
                                                                       5638, 9163, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 13738, 0, 3,
                                                                       12583, 8563, 12793, 5638,
                                                                       5788, 9388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 14053, 0, 3,
                                                                       12793, 8713, 13003, 5788,
                                                                       5938, 9613, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 14368, 0, 3,
                                                                       13003, 8863, 13213, 5938,
                                                                       6088, 9838, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 14683, 0, 3,
                                                                       13423, 9163, 13738, 6388,
                                                                       6598, 10063, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 15124, 0, 3,
                                                                       13738, 9388, 14053, 6598,
                                                                       6808, 10378, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 15565, 0, 3,
                                                                       14053, 9613, 14368, 6808,
                                                                       7018, 10693, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16006, 3, 7438,
                                                                       7453, 11050, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16034, 3, 7453,
                                                                       7468, 11071, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16062, 3, 7468,
                                                                       7483, 11092, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16090, 3, 7483,
                                                                       7498, 11113, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16118, 3, 7498,
                                                                       7513, 11134, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16146, 3, 7513,
                                                                       7528, 11155, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16174, 0, 3,
                                                                       16006, 11050, 16034, 7558,
                                                                       7603, 11302, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16258, 0, 3,
                                                                       16034, 11071, 16062, 7603,
                                                                       7648, 11365, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16342, 0, 3,
                                                                       16062, 11092, 16090, 7648,
                                                                       7693, 11428, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16426, 0, 3,
                                                                       16090, 11113, 16118, 7693,
                                                                       7738, 11491, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16510, 0, 3,
                                                                       16118, 11134, 16146, 7738,
                                                                       7783, 11554, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 16594, 0, 3,
                                                                       16174, 11302, 16258, 7873,
                                                                       7963, 11869, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 16762, 0, 3,
                                                                       16258, 11365, 16342, 7963,
                                                                       8053, 11995, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 16930, 0, 3,
                                                                       16342, 11428, 16426, 8053,
                                                                       8143, 12121, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 17098, 0, 3,
                                                                       16426, 11491, 16510, 8143,
                                                                       8233, 12247, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 17266, 0, 3,
                                                                       16594, 11869, 16762, 8413,
                                                                       8563, 12793, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 17546, 0, 3,
                                                                       16762, 11995, 16930, 8563,
                                                                       8713, 13003, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 17826, 0, 3,
                                                                       16930, 12121, 17098, 8713,
                                                                       8863, 13213, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 18106, 0, 3,
                                                                       17266, 12793, 17546, 9163,
                                                                       9388, 14053, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 18526, 0, 3,
                                                                       17546, 13003, 17826, 9388,
                                                                       9613, 14368, ncols, gamma,
                                                                       p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 18946, 0, 3,
                                                                       18106, 14053, 18526,
                                                                       10063, 10378, 15565,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19534, 3, 11008,
                                                                       11029, 16006, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19570, 3, 11029,
                                                                       11050, 16034, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19606, 3, 11050,
                                                                       11071, 16062, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19642, 3, 11071,
                                                                       11092, 16090, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19678, 3, 11092,
                                                                       11113, 16118, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19714, 3, 11113,
                                                                       11134, 16146, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 19750, 0, 3,
                                                                       19534, 16006, 19570,
                                                                       11176, 11239, 16174,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 19858, 0, 3,
                                                                       19570, 16034, 19606,
                                                                       11239, 11302, 16258,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 19966, 0, 3,
                                                                       19606, 16062, 19642,
                                                                       11302, 11365, 16342,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 20074, 0, 3,
                                                                       19642, 16090, 19678,
                                                                       11365, 11428, 16426,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 20182, 0, 3,
                                                                       19678, 16118, 19714,
                                                                       11428, 11491, 16510,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 20290, 0, 3,
                                                                       19750, 16174, 19858,
                                                                       11617, 11743, 16594,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 20506, 0, 3,
                                                                       19858, 16258, 19966,
                                                                       11743, 11869, 16762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 20722, 0, 3,
                                                                       19966, 16342, 20074,
                                                                       11869, 11995, 16930,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 20938, 0, 3,
                                                                       20074, 16426, 20182,
                                                                       11995, 12121, 17098,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 21154, 0, 3,
                                                                       20290, 16594, 20506,
                                                                       12373, 12583, 17266,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 21514, 0, 3,
                                                                       20506, 16762, 20722,
                                                                       12583, 12793, 17546,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 21874, 0, 3,
                                                                       20722, 16930, 20938,
                                                                       12793, 13003, 17826,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 22234, 0, 3,
                                                                       21154, 17266, 21514,
                                                                       13423, 13738, 18106,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 22774, 0, 3,
                                                                       21514, 17546, 21874,
                                                                       13738, 14053, 18526,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 23314, 0, 3,
                                                                       22234, 18106, 22774,
                                                                       14683, 15124, 18946,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 24070, 22234, 540, ncols);

                    simdfunc::contract_primitives(buffer, 24835, 23314, 756, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 24610, 24070, 15, 1, nmax);

        simdtrf::transform_k_inner(buffer, 25591, 24835, 21, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 25906, 24610, 25591, 15, nmax);

        simdtrf::transform_g_inner(buffer, 26581, 25906, 3, 15, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 26581, 135, nmax);
    }

    for (size_t m = 0; m < 405; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
