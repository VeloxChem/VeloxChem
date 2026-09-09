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

    const auto nmax = simdfunc::prepare_buffer(buffer, 67532, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 67532, 30201, 4690, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 6, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
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

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 581, 0, 3, 242,
                                                                       257, 392, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 609, 0, 3, 257,
                                                                       272, 413, 434, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 637, 0, 3, 272,
                                                                       287, 434, 455, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 665, 0, 3, 287,
                                                                       302, 455, 476, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 693, 0, 3, 302,
                                                                       317, 476, 497, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 721, 0, 3, 317,
                                                                       332, 497, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 749, 0, 3, 332,
                                                                       347, 518, 539, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 777, 0, 3, 347,
                                                                       362, 539, 560, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 805, 0, 3, 392,
                                                                       413, 581, 609, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 841, 0, 3, 413,
                                                                       434, 609, 637, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 877, 0, 3, 434,
                                                                       455, 637, 665, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 913, 0, 3, 455,
                                                                       476, 665, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 949, 0, 3, 476,
                                                                       497, 693, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 985, 0, 3, 497,
                                                                       518, 721, 749, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1021, 0, 3, 518,
                                                                       539, 749, 777, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1057, 0, 3, 581,
                                                                       609, 805, 841, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1102, 0, 3, 609,
                                                                       637, 841, 877, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1147, 0, 3, 637,
                                                                       665, 877, 913, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1192, 0, 3, 665,
                                                                       693, 913, 949, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1237, 0, 3, 693,
                                                                       721, 949, 985, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1282, 0, 3, 721,
                                                                       749, 985, 1021, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1327, 0, 3, 805,
                                                                       841, 1057, 1102, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1382, 0, 3, 841,
                                                                       877, 1102, 1147, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1437, 0, 3, 877,
                                                                       913, 1147, 1192, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 913,
                                                                       949, 1192, 1237, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1547, 0, 3, 949,
                                                                       985, 1237, 1282, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1602, 0, 3, 1057,
                                                                       1102, 1327, 1382, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1668, 0, 3, 1102,
                                                                       1147, 1382, 1437, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1734, 0, 3, 1147,
                                                                       1192, 1437, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 1192,
                                                                       1237, 1492, 1547, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 1866, 0, 3, 1327,
                                                                       1382, 1602, 1668, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 1944, 0, 3, 1382,
                                                                       1437, 1668, 1734, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2022, 0, 3, 1437,
                                                                       1492, 1734, 1800, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2100, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2103, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2106, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2109, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2112, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2115, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2118, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2121, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2124, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2127, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2130, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2133, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2136, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2139, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2142, 3, 9, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2151, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2160, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2169, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2178, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2187, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2196, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2205, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2214, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2223, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2232, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2241, 3, 21, 60,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2259, 3, 24, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2277, 3, 27, 72,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2295, 3, 30, 78,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2313, 3, 33, 84,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2331, 3, 36, 90,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2349, 3, 39, 96,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2367, 3, 42, 102,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2385, 3, 45, 108,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2403, 3, 48, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2421, 3, 51, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2439, 3, 54, 126,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2457, 3, 60, 132,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2487, 3, 66, 142,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2517, 3, 72, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2547, 3, 78, 162,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2577, 3, 84, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2607, 3, 90, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2637, 3, 96, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2667, 3, 102, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2697, 3, 108, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2727, 3, 114, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2757, 3, 120, 232,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2787, 3, 132, 242,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2832, 3, 142, 257,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2877, 3, 152, 272,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2922, 3, 162, 287,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2967, 3, 172, 302,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3012, 3, 182, 317,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3057, 3, 192, 332,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3102, 3, 202, 347,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3147, 3, 212, 362,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3192, 3, 222, 377,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3237, 3, 242, 392,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3300, 3, 257, 413,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3363, 3, 272, 434,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3426, 3, 287, 455,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3489, 3, 302, 476,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3552, 3, 317, 497,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3615, 3, 332, 518,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3678, 3, 347, 539,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3741, 3, 362, 560,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3804, 3, 392, 581,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3888, 3, 413, 609,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3972, 3, 434, 637,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4056, 3, 455, 665,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4140, 3, 476, 693,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4224, 3, 497, 721,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4308, 3, 518, 749,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4392, 3, 539, 777,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4476, 3, 581, 805,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4584, 3, 609, 841,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4692, 3, 637, 877,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4800, 3, 665, 913,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4908, 3, 693, 949,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5016, 3, 721, 985,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5124, 3, 749,
                                                                       1021, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5232, 3, 805,
                                                                       1057, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5367, 3, 841,
                                                                       1102, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5502, 3, 877,
                                                                       1147, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5637, 3, 913,
                                                                       1192, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5772, 3, 949,
                                                                       1237, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5907, 3, 985,
                                                                       1282, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6042, 3, 1057,
                                                                       1327, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6207, 3, 1102,
                                                                       1382, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6372, 3, 1147,
                                                                       1437, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6537, 3, 1192,
                                                                       1492, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6702, 3, 1237,
                                                                       1547, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 6867, 3, 1327,
                                                                       1602, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7065, 3, 1382,
                                                                       1668, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7263, 3, 1437,
                                                                       1734, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7461, 3, 1492,
                                                                       1800, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 7659, 3, 1602,
                                                                       1866, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 7893, 3, 1668,
                                                                       1944, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 8127, 3, 1734,
                                                                       2022, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8361, 3, 7, 8,
                                                                       2106, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8367, 3, 8, 9,
                                                                       2109, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8373, 3, 9, 10,
                                                                       2112, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8379, 3, 10, 11,
                                                                       2115, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8385, 3, 11, 12,
                                                                       2118, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8391, 3, 12, 13,
                                                                       2121, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8397, 3, 13, 14,
                                                                       2124, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8403, 3, 14, 15,
                                                                       2127, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8409, 3, 15, 16,
                                                                       2130, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8415, 3, 16, 17,
                                                                       2133, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8421, 3, 17, 18,
                                                                       2136, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8427, 3, 18, 19,
                                                                       2139, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8433, 0, 3, 8361,
                                                                       2106, 8367, 2142, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8451, 0, 3, 8367,
                                                                       2109, 8373, 2151, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8469, 0, 3, 8373,
                                                                       2112, 8379, 2160, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8487, 0, 3, 8379,
                                                                       2115, 8385, 2169, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8505, 0, 3, 8385,
                                                                       2118, 8391, 2178, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8523, 0, 3, 8391,
                                                                       2121, 8397, 2187, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8541, 0, 3, 8397,
                                                                       2124, 8403, 2196, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8559, 0, 3, 8403,
                                                                       2127, 8409, 2205, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8577, 0, 3, 8409,
                                                                       2130, 8415, 2214, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8595, 0, 3, 8415,
                                                                       2133, 8421, 2223, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8613, 0, 3, 8421,
                                                                       2136, 8427, 2232, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8631, 0, 3, 8433,
                                                                       2142, 8451, 60, 66, 2277,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8667, 0, 3, 8451,
                                                                       2151, 8469, 66, 72, 2295,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8703, 0, 3, 8469,
                                                                       2160, 8487, 72, 78, 2313,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8739, 0, 3, 8487,
                                                                       2169, 8505, 78, 84, 2331,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8775, 0, 3, 8505,
                                                                       2178, 8523, 84, 90, 2349,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8811, 0, 3, 8523,
                                                                       2187, 8541, 90, 96, 2367,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8847, 0, 3, 8541,
                                                                       2196, 8559, 96, 102, 2385,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8883, 0, 3, 8559,
                                                                       2205, 8577, 102, 108,
                                                                       2403, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8919, 0, 3, 8577,
                                                                       2214, 8595, 108, 114,
                                                                       2421, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8955, 0, 3, 8595,
                                                                       2223, 8613, 114, 120,
                                                                       2439, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8991, 0, 3, 8631,
                                                                       2277, 8667, 132, 142,
                                                                       2517, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9051, 0, 3, 8667,
                                                                       2295, 8703, 142, 152,
                                                                       2547, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9111, 0, 3, 8703,
                                                                       2313, 8739, 152, 162,
                                                                       2577, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9171, 0, 3, 8739,
                                                                       2331, 8775, 162, 172,
                                                                       2607, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9231, 0, 3, 8775,
                                                                       2349, 8811, 172, 182,
                                                                       2637, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9291, 0, 3, 8811,
                                                                       2367, 8847, 182, 192,
                                                                       2667, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9351, 0, 3, 8847,
                                                                       2385, 8883, 192, 202,
                                                                       2697, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9411, 0, 3, 8883,
                                                                       2403, 8919, 202, 212,
                                                                       2727, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9471, 0, 3, 8919,
                                                                       2421, 8955, 212, 222,
                                                                       2757, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9531, 0, 3, 8991,
                                                                       2517, 9051, 242, 257,
                                                                       2877, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9621, 0, 3, 9051,
                                                                       2547, 9111, 257, 272,
                                                                       2922, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9711, 0, 3, 9111,
                                                                       2577, 9171, 272, 287,
                                                                       2967, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9801, 0, 3, 9171,
                                                                       2607, 9231, 287, 302,
                                                                       3012, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9891, 0, 3, 9231,
                                                                       2637, 9291, 302, 317,
                                                                       3057, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9981, 0, 3, 9291,
                                                                       2667, 9351, 317, 332,
                                                                       3102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10071, 0, 3, 9351,
                                                                       2697, 9411, 332, 347,
                                                                       3147, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10161, 0, 3, 9411,
                                                                       2727, 9471, 347, 362,
                                                                       3192, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10251, 0, 3, 9531,
                                                                       2877, 9621, 392, 413,
                                                                       3363, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10377, 0, 3, 9621,
                                                                       2922, 9711, 413, 434,
                                                                       3426, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10503, 0, 3, 9711,
                                                                       2967, 9801, 434, 455,
                                                                       3489, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10629, 0, 3, 9801,
                                                                       3012, 9891, 455, 476,
                                                                       3552, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10755, 0, 3, 9891,
                                                                       3057, 9981, 476, 497,
                                                                       3615, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10881, 0, 3, 9981,
                                                                       3102, 10071, 497, 518,
                                                                       3678, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11007, 0, 3,
                                                                       10071, 3147, 10161, 518,
                                                                       539, 3741, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11133, 0, 3,
                                                                       10251, 3363, 10377, 581,
                                                                       609, 3972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11301, 0, 3,
                                                                       10377, 3426, 10503, 609,
                                                                       637, 4056, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11469, 0, 3,
                                                                       10503, 3489, 10629, 637,
                                                                       665, 4140, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11637, 0, 3,
                                                                       10629, 3552, 10755, 665,
                                                                       693, 4224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11805, 0, 3,
                                                                       10755, 3615, 10881, 693,
                                                                       721, 4308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11973, 0, 3,
                                                                       10881, 3678, 11007, 721,
                                                                       749, 4392, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12141, 0, 3,
                                                                       11133, 3972, 11301, 805,
                                                                       841, 4692, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12357, 0, 3,
                                                                       11301, 4056, 11469, 841,
                                                                       877, 4800, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12573, 0, 3,
                                                                       11469, 4140, 11637, 877,
                                                                       913, 4908, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12789, 0, 3,
                                                                       11637, 4224, 11805, 913,
                                                                       949, 5016, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13005, 0, 3,
                                                                       11805, 4308, 11973, 949,
                                                                       985, 5124, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13221, 0, 3,
                                                                       12141, 4692, 12357, 1057,
                                                                       1102, 5502, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13491, 0, 3,
                                                                       12357, 4800, 12573, 1102,
                                                                       1147, 5637, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13761, 0, 3,
                                                                       12573, 4908, 12789, 1147,
                                                                       1192, 5772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14031, 0, 3,
                                                                       12789, 5016, 13005, 1192,
                                                                       1237, 5907, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 14301, 0, 3,
                                                                       13221, 5502, 13491, 1327,
                                                                       1382, 6372, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 14631, 0, 3,
                                                                       13491, 5637, 13761, 1382,
                                                                       1437, 6537, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 14961, 0, 3,
                                                                       13761, 5772, 14031, 1437,
                                                                       1492, 6702, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 15291, 0, 3,
                                                                       14301, 6372, 14631, 1602,
                                                                       1668, 7263, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 15687, 0, 3,
                                                                       14631, 6537, 14961, 1668,
                                                                       1734, 7461, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 16083, 0, 3,
                                                                       15291, 7263, 15687, 1866,
                                                                       1944, 8127, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16551, 3, 2100,
                                                                       2103, 8361, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16561, 3, 2103,
                                                                       2106, 8367, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16571, 3, 2106,
                                                                       2109, 8373, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16581, 3, 2109,
                                                                       2112, 8379, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16591, 3, 2112,
                                                                       2115, 8385, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16601, 3, 2115,
                                                                       2118, 8391, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16611, 3, 2118,
                                                                       2121, 8397, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16621, 3, 2121,
                                                                       2124, 8403, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16631, 3, 2124,
                                                                       2127, 8409, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16641, 3, 2127,
                                                                       2130, 8415, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16651, 3, 2130,
                                                                       2133, 8421, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16661, 3, 2133,
                                                                       2136, 8427, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16671, 0, 3,
                                                                       16551, 8361, 16561, 8433,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16701, 0, 3,
                                                                       16561, 8367, 16571, 8451,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16731, 0, 3,
                                                                       16571, 8373, 16581, 8469,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16761, 0, 3,
                                                                       16581, 8379, 16591, 8487,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16791, 0, 3,
                                                                       16591, 8385, 16601, 8505,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16821, 0, 3,
                                                                       16601, 8391, 16611, 8523,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16851, 0, 3,
                                                                       16611, 8397, 16621, 8541,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16881, 0, 3,
                                                                       16621, 8403, 16631, 8559,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16911, 0, 3,
                                                                       16631, 8409, 16641, 8577,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16941, 0, 3,
                                                                       16641, 8415, 16651, 8595,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16971, 0, 3,
                                                                       16651, 8421, 16661, 8613,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17001, 0, 3,
                                                                       16671, 8433, 16701, 2241,
                                                                       2259, 8631, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17061, 0, 3,
                                                                       16701, 8451, 16731, 2259,
                                                                       2277, 8667, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17121, 0, 3,
                                                                       16731, 8469, 16761, 2277,
                                                                       2295, 8703, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17181, 0, 3,
                                                                       16761, 8487, 16791, 2295,
                                                                       2313, 8739, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17241, 0, 3,
                                                                       16791, 8505, 16821, 2313,
                                                                       2331, 8775, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17301, 0, 3,
                                                                       16821, 8523, 16851, 2331,
                                                                       2349, 8811, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17361, 0, 3,
                                                                       16851, 8541, 16881, 2349,
                                                                       2367, 8847, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17421, 0, 3,
                                                                       16881, 8559, 16911, 2367,
                                                                       2385, 8883, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17481, 0, 3,
                                                                       16911, 8577, 16941, 2385,
                                                                       2403, 8919, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17541, 0, 3,
                                                                       16941, 8595, 16971, 2403,
                                                                       2421, 8955, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17601, 0, 3,
                                                                       17001, 8631, 17061, 2457,
                                                                       2487, 8991, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17701, 0, 3,
                                                                       17061, 8667, 17121, 2487,
                                                                       2517, 9051, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17801, 0, 3,
                                                                       17121, 8703, 17181, 2517,
                                                                       2547, 9111, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17901, 0, 3,
                                                                       17181, 8739, 17241, 2547,
                                                                       2577, 9171, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18001, 0, 3,
                                                                       17241, 8775, 17301, 2577,
                                                                       2607, 9231, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18101, 0, 3,
                                                                       17301, 8811, 17361, 2607,
                                                                       2637, 9291, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18201, 0, 3,
                                                                       17361, 8847, 17421, 2637,
                                                                       2667, 9351, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18301, 0, 3,
                                                                       17421, 8883, 17481, 2667,
                                                                       2697, 9411, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18401, 0, 3,
                                                                       17481, 8919, 17541, 2697,
                                                                       2727, 9471, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18501, 0, 3,
                                                                       17601, 8991, 17701, 2787,
                                                                       2832, 9531, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18651, 0, 3,
                                                                       17701, 9051, 17801, 2832,
                                                                       2877, 9621, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18801, 0, 3,
                                                                       17801, 9111, 17901, 2877,
                                                                       2922, 9711, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18951, 0, 3,
                                                                       17901, 9171, 18001, 2922,
                                                                       2967, 9801, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19101, 0, 3,
                                                                       18001, 9231, 18101, 2967,
                                                                       3012, 9891, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19251, 0, 3,
                                                                       18101, 9291, 18201, 3012,
                                                                       3057, 9981, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19401, 0, 3,
                                                                       18201, 9351, 18301, 3057,
                                                                       3102, 10071, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19551, 0, 3,
                                                                       18301, 9411, 18401, 3102,
                                                                       3147, 10161, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19701, 0, 3,
                                                                       18501, 9531, 18651, 3237,
                                                                       3300, 10251, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19911, 0, 3,
                                                                       18651, 9621, 18801, 3300,
                                                                       3363, 10377, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20121, 0, 3,
                                                                       18801, 9711, 18951, 3363,
                                                                       3426, 10503, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20331, 0, 3,
                                                                       18951, 9801, 19101, 3426,
                                                                       3489, 10629, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20541, 0, 3,
                                                                       19101, 9891, 19251, 3489,
                                                                       3552, 10755, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20751, 0, 3,
                                                                       19251, 9981, 19401, 3552,
                                                                       3615, 10881, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20961, 0, 3,
                                                                       19401, 10071, 19551, 3615,
                                                                       3678, 11007, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21171, 0, 3,
                                                                       19701, 10251, 19911, 3804,
                                                                       3888, 11133, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21451, 0, 3,
                                                                       19911, 10377, 20121, 3888,
                                                                       3972, 11301, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21731, 0, 3,
                                                                       20121, 10503, 20331, 3972,
                                                                       4056, 11469, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22011, 0, 3,
                                                                       20331, 10629, 20541, 4056,
                                                                       4140, 11637, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22291, 0, 3,
                                                                       20541, 10755, 20751, 4140,
                                                                       4224, 11805, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22571, 0, 3,
                                                                       20751, 10881, 20961, 4224,
                                                                       4308, 11973, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 22851, 0, 3,
                                                                       21171, 11133, 21451, 4476,
                                                                       4584, 12141, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23211, 0, 3,
                                                                       21451, 11301, 21731, 4584,
                                                                       4692, 12357, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23571, 0, 3,
                                                                       21731, 11469, 22011, 4692,
                                                                       4800, 12573, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23931, 0, 3,
                                                                       22011, 11637, 22291, 4800,
                                                                       4908, 12789, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 24291, 0, 3,
                                                                       22291, 11805, 22571, 4908,
                                                                       5016, 13005, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 24651, 0, 3,
                                                                       22851, 12141, 23211, 5232,
                                                                       5367, 13221, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 25101, 0, 3,
                                                                       23211, 12357, 23571, 5367,
                                                                       5502, 13491, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 25551, 0, 3,
                                                                       23571, 12573, 23931, 5502,
                                                                       5637, 13761, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 26001, 0, 3,
                                                                       23931, 12789, 24291, 5637,
                                                                       5772, 14031, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 26451, 0, 3,
                                                                       24651, 13221, 25101, 6042,
                                                                       6207, 14301, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 27001, 0, 3,
                                                                       25101, 13491, 25551, 6207,
                                                                       6372, 14631, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 27551, 0, 3,
                                                                       25551, 13761, 26001, 6372,
                                                                       6537, 14961, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 28101, 0, 3,
                                                                       26451, 14301, 27001, 6867,
                                                                       7065, 15291, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 28761, 0, 3,
                                                                       27001, 14631, 27551, 7065,
                                                                       7263, 15687, ncols, gamma,
                                                                       p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 29421, 0, 3,
                                                                       28101, 15291, 28761, 7659,
                                                                       7893, 16083, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 30201, 21171, 280, ncols);

                    simdfunc::contract_primitives(buffer, 30677, 22851, 360, ncols);

                    simdfunc::contract_primitives(buffer, 31289, 24651, 450, ncols);

                    simdfunc::contract_primitives(buffer, 32054, 26451, 550, ncols);

                    simdfunc::contract_primitives(buffer, 32989, 28101, 660, ncols);

                    simdfunc::contract_primitives(buffer, 34111, 29421, 780, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 30481, 30201, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31037, 30677, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 31739, 31289, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 32604, 32054, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 33649, 32989, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 34891, 34111, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 35437, 30481, 31037, 7, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 36025, 31037, 31739, 7, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 36781, 31739, 32604, 7, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 37726, 32604, 33649, 7, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 38881, 33649, 34891, 7, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 40267, 35437, 36025, 7, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 41443, 36025, 36781, 7, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 42955, 36781, 37726, 7, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 44845, 37726, 38881, 7, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 47155, 40267, 41443, 7, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 49115, 41443, 42955, 7, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 51635, 42955, 44845, 7, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 54785, 47155, 49115, 7, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 57725, 49115, 51635, 7, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 61505, 54785, 57725, 7, nmax);

        simdtrf::transform_i_inner(buffer, 65621, 61505, 21, 7, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 65621, 91, nmax);
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
