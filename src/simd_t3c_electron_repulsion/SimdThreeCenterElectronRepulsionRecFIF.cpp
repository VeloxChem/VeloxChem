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


#include "SimdThreeCenterElectronRepulsionRecFIF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_fif_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_fif_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 26685, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 637 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 26685, 16050, 2403, dimensions);

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

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 818, 0, 3, 470,
                                                                       498, 638, 674, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 863, 0, 3, 498,
                                                                       526, 674, 710, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 908, 0, 3, 526,
                                                                       554, 710, 746, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 953, 0, 3, 554,
                                                                       582, 746, 782, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 998, 0, 3, 638,
                                                                       674, 818, 863, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1053, 0, 3, 674,
                                                                       710, 863, 908, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1108, 0, 3, 710,
                                                                       746, 908, 953, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1163, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1166, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1169, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1172, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1175, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1178, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1181, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1184, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1187, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1190, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1193, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1196, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1199, 3, 10, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1208, 3, 11, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1217, 3, 12, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1226, 3, 13, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1235, 3, 14, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1244, 3, 15, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1253, 3, 16, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1262, 3, 17, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1271, 3, 18, 50,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1280, 3, 20, 53,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1298, 3, 23, 59,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1316, 3, 26, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1334, 3, 29, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1352, 3, 32, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1370, 3, 35, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1388, 3, 38, 89,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1406, 3, 41, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1424, 3, 44, 101,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1442, 3, 47, 107,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1460, 3, 53, 113,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1490, 3, 59, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1520, 3, 65, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1550, 3, 71, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1580, 3, 77, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1610, 3, 83, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1640, 3, 89, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1670, 3, 95, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1700, 3, 101, 193,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1730, 3, 113, 203,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1775, 3, 123, 218,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1820, 3, 133, 233,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1865, 3, 143, 248,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1910, 3, 153, 263,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1955, 3, 163, 278,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2000, 3, 173, 293,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2045, 3, 183, 308,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2090, 3, 203, 323,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2153, 3, 218, 344,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2216, 3, 233, 365,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2279, 3, 248, 386,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2342, 3, 263, 407,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2405, 3, 278, 428,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2468, 3, 293, 449,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2531, 3, 323, 470,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2615, 3, 344, 498,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2699, 3, 365, 526,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2783, 3, 386, 554,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2867, 3, 407, 582,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2951, 3, 428, 610,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3035, 3, 470, 638,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3143, 3, 498, 674,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3251, 3, 526, 710,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3359, 3, 554, 746,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3467, 3, 582, 782,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3575, 3, 638, 818,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3710, 3, 674, 863,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3845, 3, 710, 908,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3980, 3, 746, 953,
                                                                       ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 4115, 3, 818, 998,
                                                                       ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 4280, 3, 863,
                                                                       1053, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 4445, 3, 908,
                                                                       1108, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4610, 3, 8, 9,
                                                                       1169, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4616, 3, 9, 10,
                                                                       1172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4622, 3, 10, 11,
                                                                       1175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4628, 3, 11, 12,
                                                                       1178, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4634, 3, 12, 13,
                                                                       1181, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4640, 3, 13, 14,
                                                                       1184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4646, 3, 14, 15,
                                                                       1187, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4652, 3, 15, 16,
                                                                       1190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4658, 3, 16, 17,
                                                                       1193, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4664, 3, 17, 18,
                                                                       1196, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4670, 0, 3, 4610,
                                                                       1169, 4616, 1199, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4688, 0, 3, 4616,
                                                                       1172, 4622, 1208, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4706, 0, 3, 4622,
                                                                       1175, 4628, 1217, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4724, 0, 3, 4628,
                                                                       1178, 4634, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4742, 0, 3, 4634,
                                                                       1181, 4640, 1235, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4760, 0, 3, 4640,
                                                                       1184, 4646, 1244, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4778, 0, 3, 4646,
                                                                       1187, 4652, 1253, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4796, 0, 3, 4652,
                                                                       1190, 4658, 1262, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4814, 0, 3, 4658,
                                                                       1193, 4664, 1271, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4832, 0, 3, 4670,
                                                                       1199, 4688, 53, 59, 1316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4868, 0, 3, 4688,
                                                                       1208, 4706, 59, 65, 1334,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4904, 0, 3, 4706,
                                                                       1217, 4724, 65, 71, 1352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4940, 0, 3, 4724,
                                                                       1226, 4742, 71, 77, 1370,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4976, 0, 3, 4742,
                                                                       1235, 4760, 77, 83, 1388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5012, 0, 3, 4760,
                                                                       1244, 4778, 83, 89, 1406,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5048, 0, 3, 4778,
                                                                       1253, 4796, 89, 95, 1424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5084, 0, 3, 4796,
                                                                       1262, 4814, 95, 101, 1442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5120, 0, 3, 4832,
                                                                       1316, 4868, 113, 123,
                                                                       1520, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5180, 0, 3, 4868,
                                                                       1334, 4904, 123, 133,
                                                                       1550, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5240, 0, 3, 4904,
                                                                       1352, 4940, 133, 143,
                                                                       1580, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5300, 0, 3, 4940,
                                                                       1370, 4976, 143, 153,
                                                                       1610, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5360, 0, 3, 4976,
                                                                       1388, 5012, 153, 163,
                                                                       1640, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5420, 0, 3, 5012,
                                                                       1406, 5048, 163, 173,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5480, 0, 3, 5048,
                                                                       1424, 5084, 173, 183,
                                                                       1700, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5540, 0, 3, 5120,
                                                                       1520, 5180, 203, 218,
                                                                       1820, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5630, 0, 3, 5180,
                                                                       1550, 5240, 218, 233,
                                                                       1865, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5720, 0, 3, 5240,
                                                                       1580, 5300, 233, 248,
                                                                       1910, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5810, 0, 3, 5300,
                                                                       1610, 5360, 248, 263,
                                                                       1955, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5900, 0, 3, 5360,
                                                                       1640, 5420, 263, 278,
                                                                       2000, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5990, 0, 3, 5420,
                                                                       1670, 5480, 278, 293,
                                                                       2045, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6080, 0, 3, 5540,
                                                                       1820, 5630, 323, 344,
                                                                       2216, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6206, 0, 3, 5630,
                                                                       1865, 5720, 344, 365,
                                                                       2279, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6332, 0, 3, 5720,
                                                                       1910, 5810, 365, 386,
                                                                       2342, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6458, 0, 3, 5810,
                                                                       1955, 5900, 386, 407,
                                                                       2405, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6584, 0, 3, 5900,
                                                                       2000, 5990, 407, 428,
                                                                       2468, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6710, 0, 3, 6080,
                                                                       2216, 6206, 470, 498,
                                                                       2699, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6878, 0, 3, 6206,
                                                                       2279, 6332, 498, 526,
                                                                       2783, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7046, 0, 3, 6332,
                                                                       2342, 6458, 526, 554,
                                                                       2867, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7214, 0, 3, 6458,
                                                                       2405, 6584, 554, 582,
                                                                       2951, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 7382, 0, 3, 6710,
                                                                       2699, 6878, 638, 674,
                                                                       3251, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 7598, 0, 3, 6878,
                                                                       2783, 7046, 674, 710,
                                                                       3359, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 7814, 0, 3, 7046,
                                                                       2867, 7214, 710, 746,
                                                                       3467, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 8030, 0, 3, 7382,
                                                                       3251, 7598, 818, 863,
                                                                       3845, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 8300, 0, 3, 7598,
                                                                       3359, 7814, 863, 908,
                                                                       3980, ncols, gamma, p,
                                                                       q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 8570, 0, 3, 8030,
                                                                       3845, 8300, 998, 1053,
                                                                       4445, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8900, 3, 1163,
                                                                       1166, 4610, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8910, 3, 1166,
                                                                       1169, 4616, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8920, 3, 1169,
                                                                       1172, 4622, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8930, 3, 1172,
                                                                       1175, 4628, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8940, 3, 1175,
                                                                       1178, 4634, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8950, 3, 1178,
                                                                       1181, 4640, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8960, 3, 1181,
                                                                       1184, 4646, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8970, 3, 1184,
                                                                       1187, 4652, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8980, 3, 1187,
                                                                       1190, 4658, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8990, 3, 1190,
                                                                       1193, 4664, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9000, 0, 3, 8900,
                                                                       4610, 8910, 4670, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9030, 0, 3, 8910,
                                                                       4616, 8920, 4688, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9060, 0, 3, 8920,
                                                                       4622, 8930, 4706, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9090, 0, 3, 8930,
                                                                       4628, 8940, 4724, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9120, 0, 3, 8940,
                                                                       4634, 8950, 4742, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9150, 0, 3, 8950,
                                                                       4640, 8960, 4760, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9180, 0, 3, 8960,
                                                                       4646, 8970, 4778, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9210, 0, 3, 8970,
                                                                       4652, 8980, 4796, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9240, 0, 3, 8980,
                                                                       4658, 8990, 4814, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9270, 0, 3, 9000,
                                                                       4670, 9030, 1280, 1298,
                                                                       4832, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9330, 0, 3, 9030,
                                                                       4688, 9060, 1298, 1316,
                                                                       4868, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9390, 0, 3, 9060,
                                                                       4706, 9090, 1316, 1334,
                                                                       4904, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9450, 0, 3, 9090,
                                                                       4724, 9120, 1334, 1352,
                                                                       4940, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9510, 0, 3, 9120,
                                                                       4742, 9150, 1352, 1370,
                                                                       4976, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9570, 0, 3, 9150,
                                                                       4760, 9180, 1370, 1388,
                                                                       5012, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9630, 0, 3, 9180,
                                                                       4778, 9210, 1388, 1406,
                                                                       5048, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9690, 0, 3, 9210,
                                                                       4796, 9240, 1406, 1424,
                                                                       5084, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9750, 0, 3, 9270,
                                                                       4832, 9330, 1460, 1490,
                                                                       5120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9850, 0, 3, 9330,
                                                                       4868, 9390, 1490, 1520,
                                                                       5180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9950, 0, 3, 9390,
                                                                       4904, 9450, 1520, 1550,
                                                                       5240, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10050, 0, 3, 9450,
                                                                       4940, 9510, 1550, 1580,
                                                                       5300, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10150, 0, 3, 9510,
                                                                       4976, 9570, 1580, 1610,
                                                                       5360, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10250, 0, 3, 9570,
                                                                       5012, 9630, 1610, 1640,
                                                                       5420, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10350, 0, 3, 9630,
                                                                       5048, 9690, 1640, 1670,
                                                                       5480, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10450, 0, 3, 9750,
                                                                       5120, 9850, 1730, 1775,
                                                                       5540, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10600, 0, 3, 9850,
                                                                       5180, 9950, 1775, 1820,
                                                                       5630, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10750, 0, 3, 9950,
                                                                       5240, 10050, 1820, 1865,
                                                                       5720, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10900, 0, 3,
                                                                       10050, 5300, 10150, 1865,
                                                                       1910, 5810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11050, 0, 3,
                                                                       10150, 5360, 10250, 1910,
                                                                       1955, 5900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11200, 0, 3,
                                                                       10250, 5420, 10350, 1955,
                                                                       2000, 5990, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11350, 0, 3,
                                                                       10450, 5540, 10600, 2090,
                                                                       2153, 6080, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11560, 0, 3,
                                                                       10600, 5630, 10750, 2153,
                                                                       2216, 6206, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11770, 0, 3,
                                                                       10750, 5720, 10900, 2216,
                                                                       2279, 6332, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11980, 0, 3,
                                                                       10900, 5810, 11050, 2279,
                                                                       2342, 6458, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12190, 0, 3,
                                                                       11050, 5900, 11200, 2342,
                                                                       2405, 6584, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 12400, 0, 3,
                                                                       11350, 6080, 11560, 2531,
                                                                       2615, 6710, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 12680, 0, 3,
                                                                       11560, 6206, 11770, 2615,
                                                                       2699, 6878, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 12960, 0, 3,
                                                                       11770, 6332, 11980, 2699,
                                                                       2783, 7046, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 13240, 0, 3,
                                                                       11980, 6458, 12190, 2783,
                                                                       2867, 7214, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 13520, 0, 3,
                                                                       12400, 6710, 12680, 3035,
                                                                       3143, 7382, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 13880, 0, 3,
                                                                       12680, 6878, 12960, 3143,
                                                                       3251, 7598, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 14240, 0, 3,
                                                                       12960, 7046, 13240, 3251,
                                                                       3359, 7814, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 14600, 0, 3,
                                                                       13520, 7382, 13880, 3575,
                                                                       3710, 8030, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 15050, 0, 3,
                                                                       13880, 7598, 14240, 3710,
                                                                       3845, 8300, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 15500, 0, 3,
                                                                       14600, 8030, 15050, 4115,
                                                                       4280, 8570, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 16050, 12400, 280, ncols);

                    simdfunc::contract_primitives(buffer, 16526, 13520, 360, ncols);

                    simdfunc::contract_primitives(buffer, 17138, 14600, 450, ncols);

                    simdfunc::contract_primitives(buffer, 17903, 15500, 550, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 16330, 16050, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 16886, 16526, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 17588, 17138, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 18453, 17903, 55, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 18838, 16330, 16886, 7, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 19426, 16886, 17588, 7, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 20182, 17588, 18453, 7, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 21127, 18838, 19426, 7, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 22303, 19426, 20182, 7, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 23815, 21127, 22303, 7, nmax);

        simdtrf::transform_i_inner(buffer, 25775, 23815, 10, 7, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 25775, 91, nmax);
    }

    for (size_t m = 0; m < 637; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
