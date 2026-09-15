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


#include "SimdThreeCenterElectronRepulsionRecFGG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDG.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferFG.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_fgg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_fgg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 22430, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 567 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 22430, 14198, 2076, dimensions);

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

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 470, 0, 3, 203,
                                                                       218, 323, 344, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 498, 0, 3, 218,
                                                                       233, 344, 365, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 526, 0, 3, 233,
                                                                       248, 365, 386, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 554, 0, 3, 248,
                                                                       263, 386, 407, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 582, 0, 3, 263,
                                                                       278, 407, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 610, 0, 3, 278,
                                                                       293, 428, 449, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 638, 0, 3, 323,
                                                                       344, 470, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 674, 0, 3, 344,
                                                                       365, 498, 526, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 710, 0, 3, 365,
                                                                       386, 526, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 746, 0, 3, 386,
                                                                       407, 554, 582, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 782, 0, 3, 407,
                                                                       428, 582, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 818, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 821, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 824, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 827, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 830, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 833, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 836, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 839, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 842, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 845, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 848, 3, 10, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 857, 3, 11, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 866, 3, 12, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 875, 3, 13, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 884, 3, 14, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 893, 3, 15, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 902, 3, 16, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 911, 3, 17, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 920, 3, 18, 50,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 929, 3, 26, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 947, 3, 29, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 965, 3, 32, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 983, 3, 35, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1001, 3, 38, 89,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1019, 3, 41, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1037, 3, 44, 101,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1055, 3, 47, 107,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1073, 3, 65, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1103, 3, 71, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1133, 3, 77, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1163, 3, 83, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1193, 3, 89, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1223, 3, 95, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1253, 3, 101, 193,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1283, 3, 133, 233,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1328, 3, 143, 248,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1373, 3, 153, 263,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1418, 3, 163, 278,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1463, 3, 173, 293,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1508, 3, 183, 308,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1553, 3, 233, 365,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1616, 3, 248, 386,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1679, 3, 263, 407,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1742, 3, 278, 428,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1805, 3, 293, 449,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1868, 3, 365, 526,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1952, 3, 386, 554,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2036, 3, 407, 582,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2120, 3, 428, 610,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2204, 3, 526, 710,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2312, 3, 554, 746,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2420, 3, 582, 782,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2528, 3, 8, 9,
                                                                       818, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2534, 3, 9, 10,
                                                                       821, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2540, 3, 10, 11,
                                                                       824, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2546, 3, 11, 12,
                                                                       827, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2552, 3, 12, 13,
                                                                       830, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2558, 3, 13, 14,
                                                                       833, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2564, 3, 14, 15,
                                                                       836, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2570, 3, 15, 16,
                                                                       839, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2576, 3, 16, 17,
                                                                       842, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2582, 3, 17, 18,
                                                                       845, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2588, 0, 3, 2528,
                                                                       818, 2534, 848, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2606, 0, 3, 2534,
                                                                       821, 2540, 857, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2624, 0, 3, 2540,
                                                                       824, 2546, 866, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2642, 0, 3, 2546,
                                                                       827, 2552, 875, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2660, 0, 3, 2552,
                                                                       830, 2558, 884, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2678, 0, 3, 2558,
                                                                       833, 2564, 893, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2696, 0, 3, 2564,
                                                                       836, 2570, 902, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2714, 0, 3, 2570,
                                                                       839, 2576, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2732, 0, 3, 2576,
                                                                       842, 2582, 920, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2750, 0, 3, 2588,
                                                                       848, 2606, 53, 59, 929,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2786, 0, 3, 2606,
                                                                       857, 2624, 59, 65, 947,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2822, 0, 3, 2624,
                                                                       866, 2642, 65, 71, 965,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2858, 0, 3, 2642,
                                                                       875, 2660, 71, 77, 983,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2894, 0, 3, 2660,
                                                                       884, 2678, 77, 83, 1001,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2930, 0, 3, 2678,
                                                                       893, 2696, 83, 89, 1019,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2966, 0, 3, 2696,
                                                                       902, 2714, 89, 95, 1037,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3002, 0, 3, 2714,
                                                                       911, 2732, 95, 101, 1055,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3038, 0, 3, 2750,
                                                                       929, 2786, 113, 123, 1073,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3098, 0, 3, 2786,
                                                                       947, 2822, 123, 133, 1103,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3158, 0, 3, 2822,
                                                                       965, 2858, 133, 143, 1133,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3218, 0, 3, 2858,
                                                                       983, 2894, 143, 153, 1163,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3278, 0, 3, 2894,
                                                                       1001, 2930, 153, 163,
                                                                       1193, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3338, 0, 3, 2930,
                                                                       1019, 2966, 163, 173,
                                                                       1223, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3398, 0, 3, 2966,
                                                                       1037, 3002, 173, 183,
                                                                       1253, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3458, 0, 3, 3038,
                                                                       1073, 3098, 203, 218,
                                                                       1283, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3548, 0, 3, 3098,
                                                                       1103, 3158, 218, 233,
                                                                       1328, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3638, 0, 3, 3158,
                                                                       1133, 3218, 233, 248,
                                                                       1373, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3728, 0, 3, 3218,
                                                                       1163, 3278, 248, 263,
                                                                       1418, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3818, 0, 3, 3278,
                                                                       1193, 3338, 263, 278,
                                                                       1463, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3908, 0, 3, 3338,
                                                                       1223, 3398, 278, 293,
                                                                       1508, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3998, 0, 3, 3458,
                                                                       1283, 3548, 323, 344,
                                                                       1553, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4124, 0, 3, 3548,
                                                                       1328, 3638, 344, 365,
                                                                       1616, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4250, 0, 3, 3638,
                                                                       1373, 3728, 365, 386,
                                                                       1679, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4376, 0, 3, 3728,
                                                                       1418, 3818, 386, 407,
                                                                       1742, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4502, 0, 3, 3818,
                                                                       1463, 3908, 407, 428,
                                                                       1805, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 4628, 0, 3, 3998,
                                                                       1553, 4124, 470, 498,
                                                                       1868, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 4796, 0, 3, 4124,
                                                                       1616, 4250, 498, 526,
                                                                       1952, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 4964, 0, 3, 4250,
                                                                       1679, 4376, 526, 554,
                                                                       2036, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5132, 0, 3, 4376,
                                                                       1742, 4502, 554, 582,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 5300, 0, 3, 4628,
                                                                       1868, 4796, 638, 674,
                                                                       2204, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 5516, 0, 3, 4796,
                                                                       1952, 4964, 674, 710,
                                                                       2312, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 5732, 0, 3, 4964,
                                                                       2036, 5132, 710, 746,
                                                                       2420, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5948, 3, 818, 821,
                                                                       2540, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5958, 3, 821, 824,
                                                                       2546, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5968, 3, 824, 827,
                                                                       2552, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5978, 3, 827, 830,
                                                                       2558, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5988, 3, 830, 833,
                                                                       2564, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5998, 3, 833, 836,
                                                                       2570, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6008, 3, 836, 839,
                                                                       2576, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6018, 3, 839, 842,
                                                                       2582, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6028, 0, 3, 5948,
                                                                       2540, 5958, 2624, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6058, 0, 3, 5958,
                                                                       2546, 5968, 2642, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6088, 0, 3, 5968,
                                                                       2552, 5978, 2660, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6118, 0, 3, 5978,
                                                                       2558, 5988, 2678, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6148, 0, 3, 5988,
                                                                       2564, 5998, 2696, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6178, 0, 3, 5998,
                                                                       2570, 6008, 2714, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6208, 0, 3, 6008,
                                                                       2576, 6018, 2732, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6238, 0, 3, 6028,
                                                                       2624, 6058, 929, 947,
                                                                       2822, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6298, 0, 3, 6058,
                                                                       2642, 6088, 947, 965,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6358, 0, 3, 6088,
                                                                       2660, 6118, 965, 983,
                                                                       2894, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6418, 0, 3, 6118,
                                                                       2678, 6148, 983, 1001,
                                                                       2930, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6478, 0, 3, 6148,
                                                                       2696, 6178, 1001, 1019,
                                                                       2966, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6538, 0, 3, 6178,
                                                                       2714, 6208, 1019, 1037,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6598, 0, 3, 6238,
                                                                       2822, 6298, 1073, 1103,
                                                                       3158, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6698, 0, 3, 6298,
                                                                       2858, 6358, 1103, 1133,
                                                                       3218, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6798, 0, 3, 6358,
                                                                       2894, 6418, 1133, 1163,
                                                                       3278, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6898, 0, 3, 6418,
                                                                       2930, 6478, 1163, 1193,
                                                                       3338, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 6998, 0, 3, 6478,
                                                                       2966, 6538, 1193, 1223,
                                                                       3398, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7098, 0, 3, 6598,
                                                                       3158, 6698, 1283, 1328,
                                                                       3638, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7248, 0, 3, 6698,
                                                                       3218, 6798, 1328, 1373,
                                                                       3728, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7398, 0, 3, 6798,
                                                                       3278, 6898, 1373, 1418,
                                                                       3818, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7548, 0, 3, 6898,
                                                                       3338, 6998, 1418, 1463,
                                                                       3908, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7698, 0, 3, 7098,
                                                                       3638, 7248, 1553, 1616,
                                                                       4250, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 7908, 0, 3, 7248,
                                                                       3728, 7398, 1616, 1679,
                                                                       4376, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8118, 0, 3, 7398,
                                                                       3818, 7548, 1679, 1742,
                                                                       4502, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 8328, 0, 3, 7698,
                                                                       4250, 7908, 1868, 1952,
                                                                       4964, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 8608, 0, 3, 7908,
                                                                       4376, 8118, 1952, 2036,
                                                                       5132, ncols, gamma, p,
                                                                       q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 8888, 0, 3, 8328,
                                                                       4964, 8608, 2204, 2312,
                                                                       5732, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9248, 3, 2528,
                                                                       2534, 5948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9263, 3, 2534,
                                                                       2540, 5958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9278, 3, 2540,
                                                                       2546, 5968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9293, 3, 2546,
                                                                       2552, 5978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9308, 3, 2552,
                                                                       2558, 5988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9323, 3, 2558,
                                                                       2564, 5998, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9338, 3, 2564,
                                                                       2570, 6008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9353, 3, 2570,
                                                                       2576, 6018, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9368, 0, 3, 9248,
                                                                       5948, 9263, 2588, 2606,
                                                                       6028, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9413, 0, 3, 9263,
                                                                       5958, 9278, 2606, 2624,
                                                                       6058, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9458, 0, 3, 9278,
                                                                       5968, 9293, 2624, 2642,
                                                                       6088, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9503, 0, 3, 9293,
                                                                       5978, 9308, 2642, 2660,
                                                                       6118, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9548, 0, 3, 9308,
                                                                       5988, 9323, 2660, 2678,
                                                                       6148, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9593, 0, 3, 9323,
                                                                       5998, 9338, 2678, 2696,
                                                                       6178, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9638, 0, 3, 9338,
                                                                       6008, 9353, 2696, 2714,
                                                                       6208, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9683, 0, 3, 9368,
                                                                       6028, 9413, 2750, 2786,
                                                                       6238, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9773, 0, 3, 9413,
                                                                       6058, 9458, 2786, 2822,
                                                                       6298, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9863, 0, 3, 9458,
                                                                       6088, 9503, 2822, 2858,
                                                                       6358, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9953, 0, 3, 9503,
                                                                       6118, 9548, 2858, 2894,
                                                                       6418, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 10043, 0, 3, 9548,
                                                                       6148, 9593, 2894, 2930,
                                                                       6478, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 10133, 0, 3, 9593,
                                                                       6178, 9638, 2930, 2966,
                                                                       6538, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10223, 0, 3, 9683,
                                                                       6238, 9773, 3038, 3098,
                                                                       6598, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10373, 0, 3, 9773,
                                                                       6298, 9863, 3098, 3158,
                                                                       6698, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10523, 0, 3, 9863,
                                                                       6358, 9953, 3158, 3218,
                                                                       6798, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10673, 0, 3, 9953,
                                                                       6418, 10043, 3218, 3278,
                                                                       6898, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10823, 0, 3,
                                                                       10043, 6478, 10133, 3278,
                                                                       3338, 6998, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 10973, 0, 3,
                                                                       10223, 6598, 10373, 3458,
                                                                       3548, 7098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 11198, 0, 3,
                                                                       10373, 6698, 10523, 3548,
                                                                       3638, 7248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 11423, 0, 3,
                                                                       10523, 6798, 10673, 3638,
                                                                       3728, 7398, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 11648, 0, 3,
                                                                       10673, 6898, 10823, 3728,
                                                                       3818, 7548, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 11873, 0, 3,
                                                                       10973, 7098, 11198, 3998,
                                                                       4124, 7698, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 12188, 0, 3,
                                                                       11198, 7248, 11423, 4124,
                                                                       4250, 7908, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 12503, 0, 3,
                                                                       11423, 7398, 11648, 4250,
                                                                       4376, 8118, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 12818, 0, 3,
                                                                       11873, 7698, 12188, 4628,
                                                                       4796, 8328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 13238, 0, 3,
                                                                       12188, 7908, 12503, 4796,
                                                                       4964, 8608, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 13658, 0, 3,
                                                                       12818, 8328, 13238, 5300,
                                                                       5516, 8888, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 14198, 10973, 225, ncols);

                    simdfunc::contract_primitives(buffer, 14558, 11873, 315, ncols);

                    simdfunc::contract_primitives(buffer, 15062, 12818, 420, ncols);

                    simdfunc::contract_primitives(buffer, 15734, 13658, 540, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 14423, 14198, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 14873, 14558, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 15482, 15062, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 16274, 15734, 36, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 16598, 14423, 14873, 9, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 17003, 14873, 15482, 9, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 17570, 15482, 16274, 9, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 18326, 16598, 17003, 9, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 19136, 17003, 17570, 9, nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 20270, 18326, 19136, 9, nmax);

        simdtrf::transform_g_inner(buffer, 21620, 20270, 10, 9, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 21620, 81, nmax);
    }

    for (size_t m = 0; m < 567; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
