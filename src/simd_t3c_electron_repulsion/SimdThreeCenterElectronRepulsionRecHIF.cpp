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


#include "SimdThreeCenterElectronRepulsionRecHIF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferDL.hpp"
#include "SimdTransferDM.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferFL.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferGK.hpp"
#include "SimdTransferHI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransferPN.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hif_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hif_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 67533, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1001 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 67533, 30202, 4690, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
                                                        ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 61, 0, 3, 8, 9,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 67, 0, 3, 9, 10,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 73, 0, 3, 10, 11,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 79, 0, 3, 11, 12,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 85, 0, 3, 12, 13,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 91, 0, 3, 13, 14,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 97, 0, 3, 14, 15,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 103, 0, 3, 15, 16,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 109, 0, 3, 16, 17,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 115, 0, 3, 17, 18,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 121, 0, 3, 18, 19,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 127, 0, 3, 19, 20,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 133, 0, 3, 22, 25,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 143, 0, 3, 25, 28,
                                                                       67, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 28, 31,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 31, 34,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 34, 37,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 37, 40,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 40, 43,
                                                                       97, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 43, 46,
                                                                       103, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 46, 49,
                                                                       109, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 49, 52,
                                                                       115, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 233, 0, 3, 52, 55,
                                                                       121, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 61, 67,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 67, 73,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 73, 79,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 79, 85,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 303, 0, 3, 85, 91,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 91, 97,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 333, 0, 3, 97,
                                                                       103, 193, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 103,
                                                                       109, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 363, 0, 3, 109,
                                                                       115, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 115,
                                                                       121, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 393, 0, 3, 133,
                                                                       143, 243, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 414, 0, 3, 143,
                                                                       153, 258, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 435, 0, 3, 153,
                                                                       163, 273, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 456, 0, 3, 163,
                                                                       173, 288, 303, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 477, 0, 3, 173,
                                                                       183, 303, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 183,
                                                                       193, 318, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 519, 0, 3, 193,
                                                                       203, 333, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 540, 0, 3, 203,
                                                                       213, 348, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 561, 0, 3, 213,
                                                                       223, 363, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 582, 0, 3, 243,
                                                                       258, 393, 414, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 610, 0, 3, 258,
                                                                       273, 414, 435, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 638, 0, 3, 273,
                                                                       288, 435, 456, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 666, 0, 3, 288,
                                                                       303, 456, 477, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 694, 0, 3, 303,
                                                                       318, 477, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 722, 0, 3, 318,
                                                                       333, 498, 519, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 750, 0, 3, 333,
                                                                       348, 519, 540, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 778, 0, 3, 348,
                                                                       363, 540, 561, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 806, 0, 3, 393,
                                                                       414, 582, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 842, 0, 3, 414,
                                                                       435, 610, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 878, 0, 3, 435,
                                                                       456, 638, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 914, 0, 3, 456,
                                                                       477, 666, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 950, 0, 3, 477,
                                                                       498, 694, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 986, 0, 3, 498,
                                                                       519, 722, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1022, 0, 3, 519,
                                                                       540, 750, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 582,
                                                                       610, 806, 842, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1103, 0, 3, 610,
                                                                       638, 842, 878, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1148, 0, 3, 638,
                                                                       666, 878, 914, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1193, 0, 3, 666,
                                                                       694, 914, 950, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1238, 0, 3, 694,
                                                                       722, 950, 986, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1283, 0, 3, 722,
                                                                       750, 986, 1022, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1328, 0, 3, 806,
                                                                       842, 1058, 1103, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1383, 0, 3, 842,
                                                                       878, 1103, 1148, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1438, 0, 3, 878,
                                                                       914, 1148, 1193, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1493, 0, 3, 914,
                                                                       950, 1193, 1238, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 950,
                                                                       986, 1238, 1283, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1603, 0, 3, 1058,
                                                                       1103, 1328, 1383, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1669, 0, 3, 1103,
                                                                       1148, 1383, 1438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1735, 0, 3, 1148,
                                                                       1193, 1438, 1493, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1801, 0, 3, 1193,
                                                                       1238, 1493, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 1867, 0, 3, 1328,
                                                                       1383, 1603, 1669, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 1945, 0, 3, 1383,
                                                                       1438, 1669, 1735, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2023, 0, 3, 1438,
                                                                       1493, 1735, 1801, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2101, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2104, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2107, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2110, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2113, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2116, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2119, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2122, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2125, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2128, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2131, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2134, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2137, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2140, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2143, 3, 10, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2152, 3, 11, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2161, 3, 12, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2170, 3, 13, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2179, 3, 14, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2188, 3, 15, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2197, 3, 16, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2206, 3, 17, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2215, 3, 18, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2224, 3, 19, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2233, 3, 20, 58,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2242, 3, 22, 61,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2260, 3, 25, 67,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2278, 3, 28, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2296, 3, 31, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2314, 3, 34, 85,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2332, 3, 37, 91,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2350, 3, 40, 97,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2368, 3, 43, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2386, 3, 46, 109,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2404, 3, 49, 115,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2422, 3, 52, 121,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2440, 3, 55, 127,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2458, 3, 61, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2488, 3, 67, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2518, 3, 73, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2548, 3, 79, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2578, 3, 85, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2608, 3, 91, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2638, 3, 97, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2668, 3, 103, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2698, 3, 109, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2728, 3, 115, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2758, 3, 121, 233,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2788, 3, 133, 243,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2833, 3, 143, 258,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2878, 3, 153, 273,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2923, 3, 163, 288,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2968, 3, 173, 303,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3013, 3, 183, 318,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3058, 3, 193, 333,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3103, 3, 203, 348,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3148, 3, 213, 363,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3193, 3, 223, 378,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3238, 3, 243, 393,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3301, 3, 258, 414,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3364, 3, 273, 435,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3427, 3, 288, 456,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3490, 3, 303, 477,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3553, 3, 318, 498,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3616, 3, 333, 519,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3679, 3, 348, 540,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3742, 3, 363, 561,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3805, 3, 393, 582,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3889, 3, 414, 610,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3973, 3, 435, 638,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4057, 3, 456, 666,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4141, 3, 477, 694,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4225, 3, 498, 722,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4309, 3, 519, 750,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4393, 3, 540, 778,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4477, 3, 582, 806,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4585, 3, 610, 842,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4693, 3, 638, 878,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4801, 3, 666, 914,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4909, 3, 694, 950,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5017, 3, 722, 986,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5125, 3, 750,
                                                                       1022, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5233, 3, 806,
                                                                       1058, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5368, 3, 842,
                                                                       1103, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5503, 3, 878,
                                                                       1148, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5638, 3, 914,
                                                                       1193, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5773, 3, 950,
                                                                       1238, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5908, 3, 986,
                                                                       1283, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6043, 3, 1058,
                                                                       1328, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6208, 3, 1103,
                                                                       1383, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6373, 3, 1148,
                                                                       1438, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6538, 3, 1193,
                                                                       1493, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6703, 3, 1238,
                                                                       1548, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 6868, 3, 1328,
                                                                       1603, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7066, 3, 1383,
                                                                       1669, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7264, 3, 1438,
                                                                       1735, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7462, 3, 1493,
                                                                       1801, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 7660, 3, 1603,
                                                                       1867, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 7894, 3, 1669,
                                                                       1945, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 8128, 3, 1735,
                                                                       2023, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8362, 3, 8, 9,
                                                                       2107, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8368, 3, 9, 10,
                                                                       2110, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8374, 3, 10, 11,
                                                                       2113, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8380, 3, 11, 12,
                                                                       2116, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8386, 3, 12, 13,
                                                                       2119, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8392, 3, 13, 14,
                                                                       2122, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8398, 3, 14, 15,
                                                                       2125, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8404, 3, 15, 16,
                                                                       2128, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8410, 3, 16, 17,
                                                                       2131, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8416, 3, 17, 18,
                                                                       2134, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8422, 3, 18, 19,
                                                                       2137, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8428, 3, 19, 20,
                                                                       2140, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8434, 0, 3, 8362,
                                                                       2107, 8368, 2143, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8452, 0, 3, 8368,
                                                                       2110, 8374, 2152, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8470, 0, 3, 8374,
                                                                       2113, 8380, 2161, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8488, 0, 3, 8380,
                                                                       2116, 8386, 2170, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8506, 0, 3, 8386,
                                                                       2119, 8392, 2179, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8524, 0, 3, 8392,
                                                                       2122, 8398, 2188, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8542, 0, 3, 8398,
                                                                       2125, 8404, 2197, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8560, 0, 3, 8404,
                                                                       2128, 8410, 2206, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8578, 0, 3, 8410,
                                                                       2131, 8416, 2215, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8596, 0, 3, 8416,
                                                                       2134, 8422, 2224, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8614, 0, 3, 8422,
                                                                       2137, 8428, 2233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8632, 0, 3, 8434,
                                                                       2143, 8452, 61, 67, 2278,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8668, 0, 3, 8452,
                                                                       2152, 8470, 67, 73, 2296,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8704, 0, 3, 8470,
                                                                       2161, 8488, 73, 79, 2314,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8740, 0, 3, 8488,
                                                                       2170, 8506, 79, 85, 2332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8776, 0, 3, 8506,
                                                                       2179, 8524, 85, 91, 2350,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8812, 0, 3, 8524,
                                                                       2188, 8542, 91, 97, 2368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8848, 0, 3, 8542,
                                                                       2197, 8560, 97, 103, 2386,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8884, 0, 3, 8560,
                                                                       2206, 8578, 103, 109,
                                                                       2404, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8920, 0, 3, 8578,
                                                                       2215, 8596, 109, 115,
                                                                       2422, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8956, 0, 3, 8596,
                                                                       2224, 8614, 115, 121,
                                                                       2440, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8992, 0, 3, 8632,
                                                                       2278, 8668, 133, 143,
                                                                       2518, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9052, 0, 3, 8668,
                                                                       2296, 8704, 143, 153,
                                                                       2548, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9112, 0, 3, 8704,
                                                                       2314, 8740, 153, 163,
                                                                       2578, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9172, 0, 3, 8740,
                                                                       2332, 8776, 163, 173,
                                                                       2608, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9232, 0, 3, 8776,
                                                                       2350, 8812, 173, 183,
                                                                       2638, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9292, 0, 3, 8812,
                                                                       2368, 8848, 183, 193,
                                                                       2668, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9352, 0, 3, 8848,
                                                                       2386, 8884, 193, 203,
                                                                       2698, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9412, 0, 3, 8884,
                                                                       2404, 8920, 203, 213,
                                                                       2728, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9472, 0, 3, 8920,
                                                                       2422, 8956, 213, 223,
                                                                       2758, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9532, 0, 3, 8992,
                                                                       2518, 9052, 243, 258,
                                                                       2878, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9622, 0, 3, 9052,
                                                                       2548, 9112, 258, 273,
                                                                       2923, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9712, 0, 3, 9112,
                                                                       2578, 9172, 273, 288,
                                                                       2968, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9802, 0, 3, 9172,
                                                                       2608, 9232, 288, 303,
                                                                       3013, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9892, 0, 3, 9232,
                                                                       2638, 9292, 303, 318,
                                                                       3058, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9982, 0, 3, 9292,
                                                                       2668, 9352, 318, 333,
                                                                       3103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10072, 0, 3, 9352,
                                                                       2698, 9412, 333, 348,
                                                                       3148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10162, 0, 3, 9412,
                                                                       2728, 9472, 348, 363,
                                                                       3193, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10252, 0, 3, 9532,
                                                                       2878, 9622, 393, 414,
                                                                       3364, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10378, 0, 3, 9622,
                                                                       2923, 9712, 414, 435,
                                                                       3427, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10504, 0, 3, 9712,
                                                                       2968, 9802, 435, 456,
                                                                       3490, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10630, 0, 3, 9802,
                                                                       3013, 9892, 456, 477,
                                                                       3553, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10756, 0, 3, 9892,
                                                                       3058, 9982, 477, 498,
                                                                       3616, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10882, 0, 3, 9982,
                                                                       3103, 10072, 498, 519,
                                                                       3679, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11008, 0, 3,
                                                                       10072, 3148, 10162, 519,
                                                                       540, 3742, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11134, 0, 3,
                                                                       10252, 3364, 10378, 582,
                                                                       610, 3973, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11302, 0, 3,
                                                                       10378, 3427, 10504, 610,
                                                                       638, 4057, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11470, 0, 3,
                                                                       10504, 3490, 10630, 638,
                                                                       666, 4141, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11638, 0, 3,
                                                                       10630, 3553, 10756, 666,
                                                                       694, 4225, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11806, 0, 3,
                                                                       10756, 3616, 10882, 694,
                                                                       722, 4309, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11974, 0, 3,
                                                                       10882, 3679, 11008, 722,
                                                                       750, 4393, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12142, 0, 3,
                                                                       11134, 3973, 11302, 806,
                                                                       842, 4693, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12358, 0, 3,
                                                                       11302, 4057, 11470, 842,
                                                                       878, 4801, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12574, 0, 3,
                                                                       11470, 4141, 11638, 878,
                                                                       914, 4909, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12790, 0, 3,
                                                                       11638, 4225, 11806, 914,
                                                                       950, 5017, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13006, 0, 3,
                                                                       11806, 4309, 11974, 950,
                                                                       986, 5125, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13222, 0, 3,
                                                                       12142, 4693, 12358, 1058,
                                                                       1103, 5503, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13492, 0, 3,
                                                                       12358, 4801, 12574, 1103,
                                                                       1148, 5638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13762, 0, 3,
                                                                       12574, 4909, 12790, 1148,
                                                                       1193, 5773, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14032, 0, 3,
                                                                       12790, 5017, 13006, 1193,
                                                                       1238, 5908, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 14302, 0, 3,
                                                                       13222, 5503, 13492, 1328,
                                                                       1383, 6373, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 14632, 0, 3,
                                                                       13492, 5638, 13762, 1383,
                                                                       1438, 6538, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 14962, 0, 3,
                                                                       13762, 5773, 14032, 1438,
                                                                       1493, 6703, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 15292, 0, 3,
                                                                       14302, 6373, 14632, 1603,
                                                                       1669, 7264, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 15688, 0, 3,
                                                                       14632, 6538, 14962, 1669,
                                                                       1735, 7462, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 16084, 0, 3,
                                                                       15292, 7264, 15688, 1867,
                                                                       1945, 8128, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16552, 3, 2101,
                                                                       2104, 8362, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16562, 3, 2104,
                                                                       2107, 8368, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16572, 3, 2107,
                                                                       2110, 8374, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16582, 3, 2110,
                                                                       2113, 8380, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16592, 3, 2113,
                                                                       2116, 8386, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16602, 3, 2116,
                                                                       2119, 8392, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16612, 3, 2119,
                                                                       2122, 8398, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16622, 3, 2122,
                                                                       2125, 8404, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16632, 3, 2125,
                                                                       2128, 8410, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16642, 3, 2128,
                                                                       2131, 8416, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16652, 3, 2131,
                                                                       2134, 8422, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16662, 3, 2134,
                                                                       2137, 8428, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16672, 0, 3,
                                                                       16552, 8362, 16562, 8434,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16702, 0, 3,
                                                                       16562, 8368, 16572, 8452,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16732, 0, 3,
                                                                       16572, 8374, 16582, 8470,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16762, 0, 3,
                                                                       16582, 8380, 16592, 8488,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16792, 0, 3,
                                                                       16592, 8386, 16602, 8506,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16822, 0, 3,
                                                                       16602, 8392, 16612, 8524,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16852, 0, 3,
                                                                       16612, 8398, 16622, 8542,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16882, 0, 3,
                                                                       16622, 8404, 16632, 8560,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16912, 0, 3,
                                                                       16632, 8410, 16642, 8578,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16942, 0, 3,
                                                                       16642, 8416, 16652, 8596,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16972, 0, 3,
                                                                       16652, 8422, 16662, 8614,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17002, 0, 3,
                                                                       16672, 8434, 16702, 2242,
                                                                       2260, 8632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17062, 0, 3,
                                                                       16702, 8452, 16732, 2260,
                                                                       2278, 8668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17122, 0, 3,
                                                                       16732, 8470, 16762, 2278,
                                                                       2296, 8704, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17182, 0, 3,
                                                                       16762, 8488, 16792, 2296,
                                                                       2314, 8740, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17242, 0, 3,
                                                                       16792, 8506, 16822, 2314,
                                                                       2332, 8776, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17302, 0, 3,
                                                                       16822, 8524, 16852, 2332,
                                                                       2350, 8812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17362, 0, 3,
                                                                       16852, 8542, 16882, 2350,
                                                                       2368, 8848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17422, 0, 3,
                                                                       16882, 8560, 16912, 2368,
                                                                       2386, 8884, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17482, 0, 3,
                                                                       16912, 8578, 16942, 2386,
                                                                       2404, 8920, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17542, 0, 3,
                                                                       16942, 8596, 16972, 2404,
                                                                       2422, 8956, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17602, 0, 3,
                                                                       17002, 8632, 17062, 2458,
                                                                       2488, 8992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17702, 0, 3,
                                                                       17062, 8668, 17122, 2488,
                                                                       2518, 9052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17802, 0, 3,
                                                                       17122, 8704, 17182, 2518,
                                                                       2548, 9112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17902, 0, 3,
                                                                       17182, 8740, 17242, 2548,
                                                                       2578, 9172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18002, 0, 3,
                                                                       17242, 8776, 17302, 2578,
                                                                       2608, 9232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18102, 0, 3,
                                                                       17302, 8812, 17362, 2608,
                                                                       2638, 9292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18202, 0, 3,
                                                                       17362, 8848, 17422, 2638,
                                                                       2668, 9352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18302, 0, 3,
                                                                       17422, 8884, 17482, 2668,
                                                                       2698, 9412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18402, 0, 3,
                                                                       17482, 8920, 17542, 2698,
                                                                       2728, 9472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18502, 0, 3,
                                                                       17602, 8992, 17702, 2788,
                                                                       2833, 9532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18652, 0, 3,
                                                                       17702, 9052, 17802, 2833,
                                                                       2878, 9622, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18802, 0, 3,
                                                                       17802, 9112, 17902, 2878,
                                                                       2923, 9712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18952, 0, 3,
                                                                       17902, 9172, 18002, 2923,
                                                                       2968, 9802, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19102, 0, 3,
                                                                       18002, 9232, 18102, 2968,
                                                                       3013, 9892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19252, 0, 3,
                                                                       18102, 9292, 18202, 3013,
                                                                       3058, 9982, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19402, 0, 3,
                                                                       18202, 9352, 18302, 3058,
                                                                       3103, 10072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19552, 0, 3,
                                                                       18302, 9412, 18402, 3103,
                                                                       3148, 10162, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19702, 0, 3,
                                                                       18502, 9532, 18652, 3238,
                                                                       3301, 10252, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19912, 0, 3,
                                                                       18652, 9622, 18802, 3301,
                                                                       3364, 10378, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20122, 0, 3,
                                                                       18802, 9712, 18952, 3364,
                                                                       3427, 10504, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20332, 0, 3,
                                                                       18952, 9802, 19102, 3427,
                                                                       3490, 10630, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20542, 0, 3,
                                                                       19102, 9892, 19252, 3490,
                                                                       3553, 10756, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20752, 0, 3,
                                                                       19252, 9982, 19402, 3553,
                                                                       3616, 10882, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20962, 0, 3,
                                                                       19402, 10072, 19552, 3616,
                                                                       3679, 11008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21172, 0, 3,
                                                                       19702, 10252, 19912, 3805,
                                                                       3889, 11134, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21452, 0, 3,
                                                                       19912, 10378, 20122, 3889,
                                                                       3973, 11302, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21732, 0, 3,
                                                                       20122, 10504, 20332, 3973,
                                                                       4057, 11470, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22012, 0, 3,
                                                                       20332, 10630, 20542, 4057,
                                                                       4141, 11638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22292, 0, 3,
                                                                       20542, 10756, 20752, 4141,
                                                                       4225, 11806, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22572, 0, 3,
                                                                       20752, 10882, 20962, 4225,
                                                                       4309, 11974, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 22852, 0, 3,
                                                                       21172, 11134, 21452, 4477,
                                                                       4585, 12142, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23212, 0, 3,
                                                                       21452, 11302, 21732, 4585,
                                                                       4693, 12358, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23572, 0, 3,
                                                                       21732, 11470, 22012, 4693,
                                                                       4801, 12574, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23932, 0, 3,
                                                                       22012, 11638, 22292, 4801,
                                                                       4909, 12790, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 24292, 0, 3,
                                                                       22292, 11806, 22572, 4909,
                                                                       5017, 13006, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 24652, 0, 3,
                                                                       22852, 12142, 23212, 5233,
                                                                       5368, 13222, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 25102, 0, 3,
                                                                       23212, 12358, 23572, 5368,
                                                                       5503, 13492, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 25552, 0, 3,
                                                                       23572, 12574, 23932, 5503,
                                                                       5638, 13762, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 26002, 0, 3,
                                                                       23932, 12790, 24292, 5638,
                                                                       5773, 14032, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 26452, 0, 3,
                                                                       24652, 13222, 25102, 6043,
                                                                       6208, 14302, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 27002, 0, 3,
                                                                       25102, 13492, 25552, 6208,
                                                                       6373, 14632, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 27552, 0, 3,
                                                                       25552, 13762, 26002, 6373,
                                                                       6538, 14962, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 28102, 0, 3,
                                                                       26452, 14302, 27002, 6868,
                                                                       7066, 15292, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 28762, 0, 3,
                                                                       27002, 14632, 27552, 7066,
                                                                       7264, 15688, ncols, gamma,
                                                                       p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 29422, 0, 3,
                                                                       28102, 15292, 28762, 7660,
                                                                       7894, 16084, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 30202, 21172, 280, ncols);

                    simdfunc::contract_primitives(buffer, 30678, 22852, 360, ncols);

                    simdfunc::contract_primitives(buffer, 31290, 24652, 450, ncols);

                    simdfunc::contract_primitives(buffer, 32055, 26452, 550, ncols);

                    simdfunc::contract_primitives(buffer, 32990, 28102, 660, ncols);

                    simdfunc::contract_primitives(buffer, 34112, 29422, 780, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 30482, 30202, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31038, 30678, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31740, 31290, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 32605, 32055, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 33650, 32990, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 34892, 34112, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 35438, 30482, 31038, 7, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 36026, 31038, 31740, 7, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 36782, 31740, 32605, 7, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 37727, 32605, 33650, 7, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 38882, 33650, 34892, 7, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 40268, 35438, 36026, 7, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 41444, 36026, 36782, 7, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 42956, 36782, 37727, 7, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 44846, 37727, 38882, 7, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 47156, 40268, 41444, 7, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 49116, 41444, 42956, 7, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 51636, 42956, 44846, 7, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 54786, 47156, 49116, 7, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 57726, 49116, 51636, 7, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 61506, 54786, 57726, 7, nmax);

        simdtrf::transform_i_inner(buffer, 65622, 61506, 21, 7, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 65622, 91, nmax);
    }

    for (size_t m = 0; m < 1001; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
