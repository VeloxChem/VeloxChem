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


#include "SimdThreeCenterElectronRepulsionRecDIF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_dif_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_dif_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 16129, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 455 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 16129, 11210, 1538, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11}, ncols, fj, 6,
                                                        fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 19, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 49, 0, 3, 8, 9,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 55, 0, 3, 9, 10,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 61, 0, 3, 10, 11,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 67, 0, 3, 11, 12,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 73, 0, 3, 12, 13,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 79, 0, 3, 13, 14,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 85, 0, 3, 14, 15,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 91, 0, 3, 15, 16,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 97, 0, 3, 16, 17,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 103, 0, 3, 19, 22,
                                                                       49, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 113, 0, 3, 22, 25,
                                                                       55, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 123, 0, 3, 25, 28,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 133, 0, 3, 28, 31,
                                                                       67, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 143, 0, 3, 31, 34,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 34, 37,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 37, 40,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 40, 43,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 49, 55,
                                                                       103, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 198, 0, 3, 55, 61,
                                                                       113, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 61, 67,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 67, 73,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 73, 79,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 79, 85,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 85, 91,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 103,
                                                                       113, 183, 198, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 309, 0, 3, 113,
                                                                       123, 198, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 330, 0, 3, 123,
                                                                       133, 213, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 351, 0, 3, 133,
                                                                       143, 228, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 372, 0, 3, 143,
                                                                       153, 243, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 393, 0, 3, 153,
                                                                       163, 258, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 414, 0, 3, 183,
                                                                       198, 288, 309, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 442, 0, 3, 198,
                                                                       213, 309, 330, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 470, 0, 3, 213,
                                                                       228, 330, 351, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 498, 0, 3, 228,
                                                                       243, 351, 372, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 526, 0, 3, 243,
                                                                       258, 372, 393, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 554, 0, 3, 288,
                                                                       309, 414, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 590, 0, 3, 309,
                                                                       330, 442, 470, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 626, 0, 3, 330,
                                                                       351, 470, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 662, 0, 3, 351,
                                                                       372, 498, 526, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 698, 0, 3, 414,
                                                                       442, 554, 590, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 743, 0, 3, 442,
                                                                       470, 590, 626, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 788, 0, 3, 470,
                                                                       498, 626, 662, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 833, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 836, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 839, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 842, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 845, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 848, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 851, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 854, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 857, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 860, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 863, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 866, 3, 10, 25,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 875, 3, 11, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 884, 3, 12, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 893, 3, 13, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 902, 3, 14, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 911, 3, 15, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 920, 3, 16, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 929, 3, 17, 46,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 938, 3, 19, 49,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 956, 3, 22, 55,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 974, 3, 25, 61,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 992, 3, 28, 67,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1010, 3, 31, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1028, 3, 34, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1046, 3, 37, 85,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1064, 3, 40, 91,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1082, 3, 43, 97,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1100, 3, 49, 103,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1130, 3, 55, 113,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1160, 3, 61, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1190, 3, 67, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1220, 3, 73, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1250, 3, 79, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1280, 3, 85, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1310, 3, 91, 173,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1340, 3, 103, 183,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1385, 3, 113, 198,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1430, 3, 123, 213,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1475, 3, 133, 228,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1520, 3, 143, 243,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1565, 3, 153, 258,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1610, 3, 163, 273,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1655, 3, 183, 288,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1718, 3, 198, 309,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1781, 3, 213, 330,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1844, 3, 228, 351,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1907, 3, 243, 372,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1970, 3, 258, 393,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2033, 3, 288, 414,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2117, 3, 309, 442,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2201, 3, 330, 470,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2285, 3, 351, 498,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2369, 3, 372, 526,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2453, 3, 414, 554,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2561, 3, 442, 590,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2669, 3, 470, 626,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2777, 3, 498, 662,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 2885, 3, 554, 698,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3020, 3, 590, 743,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3155, 3, 626, 788,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3290, 3, 8, 9,
                                                                       839, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3296, 3, 9, 10,
                                                                       842, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3302, 3, 10, 11,
                                                                       845, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3308, 3, 11, 12,
                                                                       848, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3314, 3, 12, 13,
                                                                       851, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3320, 3, 13, 14,
                                                                       854, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3326, 3, 14, 15,
                                                                       857, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3332, 3, 15, 16,
                                                                       860, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3338, 3, 16, 17,
                                                                       863, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3344, 0, 3, 3290,
                                                                       839, 3296, 866, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3362, 0, 3, 3296,
                                                                       842, 3302, 875, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3380, 0, 3, 3302,
                                                                       845, 3308, 884, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3398, 0, 3, 3308,
                                                                       848, 3314, 893, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3416, 0, 3, 3314,
                                                                       851, 3320, 902, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3434, 0, 3, 3320,
                                                                       854, 3326, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3452, 0, 3, 3326,
                                                                       857, 3332, 920, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3470, 0, 3, 3332,
                                                                       860, 3338, 929, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3488, 0, 3, 3344,
                                                                       866, 3362, 49, 55, 974,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3524, 0, 3, 3362,
                                                                       875, 3380, 55, 61, 992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3560, 0, 3, 3380,
                                                                       884, 3398, 61, 67, 1010,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3596, 0, 3, 3398,
                                                                       893, 3416, 67, 73, 1028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3632, 0, 3, 3416,
                                                                       902, 3434, 73, 79, 1046,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3668, 0, 3, 3434,
                                                                       911, 3452, 79, 85, 1064,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3704, 0, 3, 3452,
                                                                       920, 3470, 85, 91, 1082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3740, 0, 3, 3488,
                                                                       974, 3524, 103, 113, 1160,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3800, 0, 3, 3524,
                                                                       992, 3560, 113, 123, 1190,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3860, 0, 3, 3560,
                                                                       1010, 3596, 123, 133,
                                                                       1220, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3920, 0, 3, 3596,
                                                                       1028, 3632, 133, 143,
                                                                       1250, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3980, 0, 3, 3632,
                                                                       1046, 3668, 143, 153,
                                                                       1280, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4040, 0, 3, 3668,
                                                                       1064, 3704, 153, 163,
                                                                       1310, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4100, 0, 3, 3740,
                                                                       1160, 3800, 183, 198,
                                                                       1430, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4190, 0, 3, 3800,
                                                                       1190, 3860, 198, 213,
                                                                       1475, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4280, 0, 3, 3860,
                                                                       1220, 3920, 213, 228,
                                                                       1520, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4370, 0, 3, 3920,
                                                                       1250, 3980, 228, 243,
                                                                       1565, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4460, 0, 3, 3980,
                                                                       1280, 4040, 243, 258,
                                                                       1610, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4550, 0, 3, 4100,
                                                                       1430, 4190, 288, 309,
                                                                       1781, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4676, 0, 3, 4190,
                                                                       1475, 4280, 309, 330,
                                                                       1844, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4802, 0, 3, 4280,
                                                                       1520, 4370, 330, 351,
                                                                       1907, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4928, 0, 3, 4370,
                                                                       1565, 4460, 351, 372,
                                                                       1970, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5054, 0, 3, 4550,
                                                                       1781, 4676, 414, 442,
                                                                       2201, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5222, 0, 3, 4676,
                                                                       1844, 4802, 442, 470,
                                                                       2285, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5390, 0, 3, 4802,
                                                                       1907, 4928, 470, 498,
                                                                       2369, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 5558, 0, 3, 5054,
                                                                       2201, 5222, 554, 590,
                                                                       2669, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 5774, 0, 3, 5222,
                                                                       2285, 5390, 590, 626,
                                                                       2777, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 5990, 0, 3, 5558,
                                                                       2669, 5774, 698, 743,
                                                                       3155, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6260, 3, 833, 836,
                                                                       3290, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6270, 3, 836, 839,
                                                                       3296, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6280, 3, 839, 842,
                                                                       3302, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6290, 3, 842, 845,
                                                                       3308, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6300, 3, 845, 848,
                                                                       3314, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6310, 3, 848, 851,
                                                                       3320, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6320, 3, 851, 854,
                                                                       3326, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6330, 3, 854, 857,
                                                                       3332, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6340, 3, 857, 860,
                                                                       3338, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6350, 0, 3, 6260,
                                                                       3290, 6270, 3344, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6380, 0, 3, 6270,
                                                                       3296, 6280, 3362, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6410, 0, 3, 6280,
                                                                       3302, 6290, 3380, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6440, 0, 3, 6290,
                                                                       3308, 6300, 3398, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6470, 0, 3, 6300,
                                                                       3314, 6310, 3416, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6500, 0, 3, 6310,
                                                                       3320, 6320, 3434, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6530, 0, 3, 6320,
                                                                       3326, 6330, 3452, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6560, 0, 3, 6330,
                                                                       3332, 6340, 3470, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6590, 0, 3, 6350,
                                                                       3344, 6380, 938, 956,
                                                                       3488, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6650, 0, 3, 6380,
                                                                       3362, 6410, 956, 974,
                                                                       3524, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6710, 0, 3, 6410,
                                                                       3380, 6440, 974, 992,
                                                                       3560, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6770, 0, 3, 6440,
                                                                       3398, 6470, 992, 1010,
                                                                       3596, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6830, 0, 3, 6470,
                                                                       3416, 6500, 1010, 1028,
                                                                       3632, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6890, 0, 3, 6500,
                                                                       3434, 6530, 1028, 1046,
                                                                       3668, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6950, 0, 3, 6530,
                                                                       3452, 6560, 1046, 1064,
                                                                       3704, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7010, 0, 3, 6590,
                                                                       3488, 6650, 1100, 1130,
                                                                       3740, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7110, 0, 3, 6650,
                                                                       3524, 6710, 1130, 1160,
                                                                       3800, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7210, 0, 3, 6710,
                                                                       3560, 6770, 1160, 1190,
                                                                       3860, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7310, 0, 3, 6770,
                                                                       3596, 6830, 1190, 1220,
                                                                       3920, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7410, 0, 3, 6830,
                                                                       3632, 6890, 1220, 1250,
                                                                       3980, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7510, 0, 3, 6890,
                                                                       3668, 6950, 1250, 1280,
                                                                       4040, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7610, 0, 3, 7010,
                                                                       3740, 7110, 1340, 1385,
                                                                       4100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7760, 0, 3, 7110,
                                                                       3800, 7210, 1385, 1430,
                                                                       4190, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7910, 0, 3, 7210,
                                                                       3860, 7310, 1430, 1475,
                                                                       4280, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8060, 0, 3, 7310,
                                                                       3920, 7410, 1475, 1520,
                                                                       4370, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8210, 0, 3, 7410,
                                                                       3980, 7510, 1520, 1565,
                                                                       4460, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8360, 0, 3, 7610,
                                                                       4100, 7760, 1655, 1718,
                                                                       4550, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8570, 0, 3, 7760,
                                                                       4190, 7910, 1718, 1781,
                                                                       4676, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8780, 0, 3, 7910,
                                                                       4280, 8060, 1781, 1844,
                                                                       4802, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8990, 0, 3, 8060,
                                                                       4370, 8210, 1844, 1907,
                                                                       4928, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 9200, 0, 3, 8360,
                                                                       4550, 8570, 2033, 2117,
                                                                       5054, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 9480, 0, 3, 8570,
                                                                       4676, 8780, 2117, 2201,
                                                                       5222, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 9760, 0, 3, 8780,
                                                                       4802, 8990, 2201, 2285,
                                                                       5390, ncols, gamma, p,
                                                                       q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 10040, 0, 3, 9200,
                                                                       5054, 9480, 2453, 2561,
                                                                       5558, ncols, gamma, p,
                                                                       q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 10400, 0, 3, 9480,
                                                                       5222, 9760, 2561, 2669,
                                                                       5774, ncols, gamma, p,
                                                                       q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 10760, 0, 3,
                                                                       10040, 5558, 10400, 2885,
                                                                       3020, 5990, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 11210, 9200, 280, ncols);

                    simdfunc::contract_primitives(buffer, 11686, 10040, 360, ncols);

                    simdfunc::contract_primitives(buffer, 12298, 10760, 450, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 11490, 11210, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 12046, 11686, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 12748, 12298, 45, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 13063, 11490, 12046, 7, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 13651, 12046, 12748, 7, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 14407, 13063, 13651, 7, nmax);

        simdtrf::transform_i_inner(buffer, 15583, 14407, 6, 7, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 15583, 91, nmax);
    }

    for (size_t m = 0; m < 455; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
