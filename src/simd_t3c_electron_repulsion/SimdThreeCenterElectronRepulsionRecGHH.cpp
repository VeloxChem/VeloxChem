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


#include "SimdThreeCenterElectronRepulsionRecGHH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferFH.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferGH.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ghh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ghh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 77079, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1089 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 77079, 50589, 5315, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1602, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1605, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1608, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1611, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1614, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1617, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1620, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1623, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1626, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1629, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1632, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1635, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1638, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1641, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1644, 3, 9, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1653, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1662, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1671, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1680, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1689, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1698, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1707, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1716, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1725, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1734, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1743, 3, 21, 60,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1761, 3, 24, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1779, 3, 27, 72,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1797, 3, 30, 78,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1815, 3, 33, 84,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1833, 3, 36, 90,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1851, 3, 39, 96,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1869, 3, 42, 102,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1887, 3, 45, 108,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1905, 3, 48, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1923, 3, 51, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1941, 3, 54, 126,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1959, 3, 60, 132,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1989, 3, 66, 142,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2019, 3, 72, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2049, 3, 78, 162,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2079, 3, 84, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2109, 3, 90, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2139, 3, 96, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2169, 3, 102, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2199, 3, 108, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2229, 3, 114, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2259, 3, 120, 232,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2289, 3, 132, 242,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2334, 3, 142, 257,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2379, 3, 152, 272,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2424, 3, 162, 287,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2469, 3, 172, 302,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2514, 3, 182, 317,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2559, 3, 192, 332,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2604, 3, 202, 347,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2649, 3, 212, 362,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2694, 3, 222, 377,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2739, 3, 242, 392,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2802, 3, 257, 413,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2865, 3, 272, 434,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2928, 3, 287, 455,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2991, 3, 302, 476,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3054, 3, 317, 497,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3117, 3, 332, 518,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3180, 3, 347, 539,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3243, 3, 362, 560,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3306, 3, 392, 581,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3390, 3, 413, 609,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3474, 3, 434, 637,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3558, 3, 455, 665,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3642, 3, 476, 693,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3726, 3, 497, 721,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3810, 3, 518, 749,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3894, 3, 539, 777,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3978, 3, 581, 805,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4086, 3, 609, 841,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4194, 3, 637, 877,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4302, 3, 665, 913,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4410, 3, 693, 949,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4518, 3, 721, 985,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4626, 3, 749,
                                                                       1021, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4734, 3, 805,
                                                                       1057, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4869, 3, 841,
                                                                       1102, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5004, 3, 877,
                                                                       1147, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5139, 3, 913,
                                                                       1192, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5274, 3, 949,
                                                                       1237, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5409, 3, 985,
                                                                       1282, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 5544, 3, 1057,
                                                                       1327, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 5709, 3, 1102,
                                                                       1382, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 5874, 3, 1147,
                                                                       1437, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6039, 3, 1192,
                                                                       1492, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6204, 3, 1237,
                                                                       1547, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6369, 3, 7, 8,
                                                                       1608, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6375, 3, 8, 9,
                                                                       1611, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6381, 3, 9, 10,
                                                                       1614, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6387, 3, 10, 11,
                                                                       1617, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6393, 3, 11, 12,
                                                                       1620, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6399, 3, 12, 13,
                                                                       1623, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6405, 3, 13, 14,
                                                                       1626, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6411, 3, 14, 15,
                                                                       1629, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6417, 3, 15, 16,
                                                                       1632, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6423, 3, 16, 17,
                                                                       1635, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6429, 3, 17, 18,
                                                                       1638, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6435, 3, 18, 19,
                                                                       1641, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6441, 0, 3, 6369,
                                                                       1608, 6375, 1644, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6459, 0, 3, 6375,
                                                                       1611, 6381, 1653, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6477, 0, 3, 6381,
                                                                       1614, 6387, 1662, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6495, 0, 3, 6387,
                                                                       1617, 6393, 1671, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6513, 0, 3, 6393,
                                                                       1620, 6399, 1680, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6531, 0, 3, 6399,
                                                                       1623, 6405, 1689, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6549, 0, 3, 6405,
                                                                       1626, 6411, 1698, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6567, 0, 3, 6411,
                                                                       1629, 6417, 1707, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6585, 0, 3, 6417,
                                                                       1632, 6423, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6603, 0, 3, 6423,
                                                                       1635, 6429, 1725, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6621, 0, 3, 6429,
                                                                       1638, 6435, 1734, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6639, 0, 3, 6441,
                                                                       1644, 6459, 60, 66, 1779,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6675, 0, 3, 6459,
                                                                       1653, 6477, 66, 72, 1797,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6711, 0, 3, 6477,
                                                                       1662, 6495, 72, 78, 1815,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6747, 0, 3, 6495,
                                                                       1671, 6513, 78, 84, 1833,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6783, 0, 3, 6513,
                                                                       1680, 6531, 84, 90, 1851,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6819, 0, 3, 6531,
                                                                       1689, 6549, 90, 96, 1869,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6855, 0, 3, 6549,
                                                                       1698, 6567, 96, 102, 1887,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6891, 0, 3, 6567,
                                                                       1707, 6585, 102, 108,
                                                                       1905, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6927, 0, 3, 6585,
                                                                       1716, 6603, 108, 114,
                                                                       1923, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6963, 0, 3, 6603,
                                                                       1725, 6621, 114, 120,
                                                                       1941, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6999, 0, 3, 6639,
                                                                       1779, 6675, 132, 142,
                                                                       2019, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7059, 0, 3, 6675,
                                                                       1797, 6711, 142, 152,
                                                                       2049, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7119, 0, 3, 6711,
                                                                       1815, 6747, 152, 162,
                                                                       2079, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7179, 0, 3, 6747,
                                                                       1833, 6783, 162, 172,
                                                                       2109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7239, 0, 3, 6783,
                                                                       1851, 6819, 172, 182,
                                                                       2139, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7299, 0, 3, 6819,
                                                                       1869, 6855, 182, 192,
                                                                       2169, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7359, 0, 3, 6855,
                                                                       1887, 6891, 192, 202,
                                                                       2199, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7419, 0, 3, 6891,
                                                                       1905, 6927, 202, 212,
                                                                       2229, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7479, 0, 3, 6927,
                                                                       1923, 6963, 212, 222,
                                                                       2259, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7539, 0, 3, 6999,
                                                                       2019, 7059, 242, 257,
                                                                       2379, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7629, 0, 3, 7059,
                                                                       2049, 7119, 257, 272,
                                                                       2424, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7719, 0, 3, 7119,
                                                                       2079, 7179, 272, 287,
                                                                       2469, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7809, 0, 3, 7179,
                                                                       2109, 7239, 287, 302,
                                                                       2514, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7899, 0, 3, 7239,
                                                                       2139, 7299, 302, 317,
                                                                       2559, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7989, 0, 3, 7299,
                                                                       2169, 7359, 317, 332,
                                                                       2604, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8079, 0, 3, 7359,
                                                                       2199, 7419, 332, 347,
                                                                       2649, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8169, 0, 3, 7419,
                                                                       2229, 7479, 347, 362,
                                                                       2694, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8259, 0, 3, 7539,
                                                                       2379, 7629, 392, 413,
                                                                       2865, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8385, 0, 3, 7629,
                                                                       2424, 7719, 413, 434,
                                                                       2928, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8511, 0, 3, 7719,
                                                                       2469, 7809, 434, 455,
                                                                       2991, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8637, 0, 3, 7809,
                                                                       2514, 7899, 455, 476,
                                                                       3054, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8763, 0, 3, 7899,
                                                                       2559, 7989, 476, 497,
                                                                       3117, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8889, 0, 3, 7989,
                                                                       2604, 8079, 497, 518,
                                                                       3180, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9015, 0, 3, 8079,
                                                                       2649, 8169, 518, 539,
                                                                       3243, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9141, 0, 3, 8259,
                                                                       2865, 8385, 581, 609,
                                                                       3474, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9309, 0, 3, 8385,
                                                                       2928, 8511, 609, 637,
                                                                       3558, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9477, 0, 3, 8511,
                                                                       2991, 8637, 637, 665,
                                                                       3642, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9645, 0, 3, 8637,
                                                                       3054, 8763, 665, 693,
                                                                       3726, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9813, 0, 3, 8763,
                                                                       3117, 8889, 693, 721,
                                                                       3810, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9981, 0, 3, 8889,
                                                                       3180, 9015, 721, 749,
                                                                       3894, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10149, 0, 3, 9141,
                                                                       3474, 9309, 805, 841,
                                                                       4194, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10365, 0, 3, 9309,
                                                                       3558, 9477, 841, 877,
                                                                       4302, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10581, 0, 3, 9477,
                                                                       3642, 9645, 877, 913,
                                                                       4410, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10797, 0, 3, 9645,
                                                                       3726, 9813, 913, 949,
                                                                       4518, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11013, 0, 3, 9813,
                                                                       3810, 9981, 949, 985,
                                                                       4626, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 11229, 0, 3,
                                                                       10149, 4194, 10365, 1057,
                                                                       1102, 5004, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 11499, 0, 3,
                                                                       10365, 4302, 10581, 1102,
                                                                       1147, 5139, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 11769, 0, 3,
                                                                       10581, 4410, 10797, 1147,
                                                                       1192, 5274, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 12039, 0, 3,
                                                                       10797, 4518, 11013, 1192,
                                                                       1237, 5409, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 12309, 0, 3,
                                                                       11229, 5004, 11499, 1327,
                                                                       1382, 5874, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 12639, 0, 3,
                                                                       11499, 5139, 11769, 1382,
                                                                       1437, 6039, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 12969, 0, 3,
                                                                       11769, 5274, 12039, 1437,
                                                                       1492, 6204, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13299, 3, 1602,
                                                                       1605, 6369, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13309, 3, 1605,
                                                                       1608, 6375, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13319, 3, 1608,
                                                                       1611, 6381, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13329, 3, 1611,
                                                                       1614, 6387, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13339, 3, 1614,
                                                                       1617, 6393, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13349, 3, 1617,
                                                                       1620, 6399, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13359, 3, 1620,
                                                                       1623, 6405, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13369, 3, 1623,
                                                                       1626, 6411, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13379, 3, 1626,
                                                                       1629, 6417, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13389, 3, 1629,
                                                                       1632, 6423, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13399, 3, 1632,
                                                                       1635, 6429, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13409, 3, 1635,
                                                                       1638, 6435, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13419, 0, 3,
                                                                       13299, 6369, 13309, 6441,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13449, 0, 3,
                                                                       13309, 6375, 13319, 6459,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13479, 0, 3,
                                                                       13319, 6381, 13329, 6477,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13509, 0, 3,
                                                                       13329, 6387, 13339, 6495,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13539, 0, 3,
                                                                       13339, 6393, 13349, 6513,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13569, 0, 3,
                                                                       13349, 6399, 13359, 6531,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13599, 0, 3,
                                                                       13359, 6405, 13369, 6549,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13629, 0, 3,
                                                                       13369, 6411, 13379, 6567,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13659, 0, 3,
                                                                       13379, 6417, 13389, 6585,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13689, 0, 3,
                                                                       13389, 6423, 13399, 6603,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13719, 0, 3,
                                                                       13399, 6429, 13409, 6621,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13749, 0, 3,
                                                                       13419, 6441, 13449, 1743,
                                                                       1761, 6639, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13809, 0, 3,
                                                                       13449, 6459, 13479, 1761,
                                                                       1779, 6675, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13869, 0, 3,
                                                                       13479, 6477, 13509, 1779,
                                                                       1797, 6711, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13929, 0, 3,
                                                                       13509, 6495, 13539, 1797,
                                                                       1815, 6747, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13989, 0, 3,
                                                                       13539, 6513, 13569, 1815,
                                                                       1833, 6783, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14049, 0, 3,
                                                                       13569, 6531, 13599, 1833,
                                                                       1851, 6819, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14109, 0, 3,
                                                                       13599, 6549, 13629, 1851,
                                                                       1869, 6855, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14169, 0, 3,
                                                                       13629, 6567, 13659, 1869,
                                                                       1887, 6891, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14229, 0, 3,
                                                                       13659, 6585, 13689, 1887,
                                                                       1905, 6927, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14289, 0, 3,
                                                                       13689, 6603, 13719, 1905,
                                                                       1923, 6963, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14349, 0, 3,
                                                                       13749, 6639, 13809, 1959,
                                                                       1989, 6999, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14449, 0, 3,
                                                                       13809, 6675, 13869, 1989,
                                                                       2019, 7059, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14549, 0, 3,
                                                                       13869, 6711, 13929, 2019,
                                                                       2049, 7119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14649, 0, 3,
                                                                       13929, 6747, 13989, 2049,
                                                                       2079, 7179, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14749, 0, 3,
                                                                       13989, 6783, 14049, 2079,
                                                                       2109, 7239, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14849, 0, 3,
                                                                       14049, 6819, 14109, 2109,
                                                                       2139, 7299, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14949, 0, 3,
                                                                       14109, 6855, 14169, 2139,
                                                                       2169, 7359, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15049, 0, 3,
                                                                       14169, 6891, 14229, 2169,
                                                                       2199, 7419, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15149, 0, 3,
                                                                       14229, 6927, 14289, 2199,
                                                                       2229, 7479, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15249, 0, 3,
                                                                       14349, 6999, 14449, 2289,
                                                                       2334, 7539, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15399, 0, 3,
                                                                       14449, 7059, 14549, 2334,
                                                                       2379, 7629, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15549, 0, 3,
                                                                       14549, 7119, 14649, 2379,
                                                                       2424, 7719, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15699, 0, 3,
                                                                       14649, 7179, 14749, 2424,
                                                                       2469, 7809, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15849, 0, 3,
                                                                       14749, 7239, 14849, 2469,
                                                                       2514, 7899, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15999, 0, 3,
                                                                       14849, 7299, 14949, 2514,
                                                                       2559, 7989, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16149, 0, 3,
                                                                       14949, 7359, 15049, 2559,
                                                                       2604, 8079, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16299, 0, 3,
                                                                       15049, 7419, 15149, 2604,
                                                                       2649, 8169, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16449, 0, 3,
                                                                       15249, 7539, 15399, 2739,
                                                                       2802, 8259, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16659, 0, 3,
                                                                       15399, 7629, 15549, 2802,
                                                                       2865, 8385, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16869, 0, 3,
                                                                       15549, 7719, 15699, 2865,
                                                                       2928, 8511, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17079, 0, 3,
                                                                       15699, 7809, 15849, 2928,
                                                                       2991, 8637, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17289, 0, 3,
                                                                       15849, 7899, 15999, 2991,
                                                                       3054, 8763, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17499, 0, 3,
                                                                       15999, 7989, 16149, 3054,
                                                                       3117, 8889, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17709, 0, 3,
                                                                       16149, 8079, 16299, 3117,
                                                                       3180, 9015, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 17919, 0, 3,
                                                                       16449, 8259, 16659, 3306,
                                                                       3390, 9141, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18199, 0, 3,
                                                                       16659, 8385, 16869, 3390,
                                                                       3474, 9309, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18479, 0, 3,
                                                                       16869, 8511, 17079, 3474,
                                                                       3558, 9477, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18759, 0, 3,
                                                                       17079, 8637, 17289, 3558,
                                                                       3642, 9645, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 19039, 0, 3,
                                                                       17289, 8763, 17499, 3642,
                                                                       3726, 9813, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 19319, 0, 3,
                                                                       17499, 8889, 17709, 3726,
                                                                       3810, 9981, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 19599, 0, 3,
                                                                       17919, 9141, 18199, 3978,
                                                                       4086, 10149, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 19959, 0, 3,
                                                                       18199, 9309, 18479, 4086,
                                                                       4194, 10365, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 20319, 0, 3,
                                                                       18479, 9477, 18759, 4194,
                                                                       4302, 10581, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 20679, 0, 3,
                                                                       18759, 9645, 19039, 4302,
                                                                       4410, 10797, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 21039, 0, 3,
                                                                       19039, 9813, 19319, 4410,
                                                                       4518, 11013, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 21399, 0, 3,
                                                                       19599, 10149, 19959, 4734,
                                                                       4869, 11229, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 21849, 0, 3,
                                                                       19959, 10365, 20319, 4869,
                                                                       5004, 11499, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 22299, 0, 3,
                                                                       20319, 10581, 20679, 5004,
                                                                       5139, 11769, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 22749, 0, 3,
                                                                       20679, 10797, 21039, 5139,
                                                                       5274, 12039, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 23199, 0, 3,
                                                                       21399, 11229, 21849, 5544,
                                                                       5709, 12309, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 23749, 0, 3,
                                                                       21849, 11499, 22299, 5709,
                                                                       5874, 12639, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 24299, 0, 3,
                                                                       22299, 11769, 22749, 5874,
                                                                       6039, 12969, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24849, 3, 6369,
                                                                       6375, 13319, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24864, 3, 6375,
                                                                       6381, 13329, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24879, 3, 6381,
                                                                       6387, 13339, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24894, 3, 6387,
                                                                       6393, 13349, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24909, 3, 6393,
                                                                       6399, 13359, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24924, 3, 6399,
                                                                       6405, 13369, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24939, 3, 6405,
                                                                       6411, 13379, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24954, 3, 6411,
                                                                       6417, 13389, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24969, 3, 6417,
                                                                       6423, 13399, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24984, 3, 6423,
                                                                       6429, 13409, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24999, 0, 3,
                                                                       24849, 13319, 24864, 6441,
                                                                       6459, 13479, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25044, 0, 3,
                                                                       24864, 13329, 24879, 6459,
                                                                       6477, 13509, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25089, 0, 3,
                                                                       24879, 13339, 24894, 6477,
                                                                       6495, 13539, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25134, 0, 3,
                                                                       24894, 13349, 24909, 6495,
                                                                       6513, 13569, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25179, 0, 3,
                                                                       24909, 13359, 24924, 6513,
                                                                       6531, 13599, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25224, 0, 3,
                                                                       24924, 13369, 24939, 6531,
                                                                       6549, 13629, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25269, 0, 3,
                                                                       24939, 13379, 24954, 6549,
                                                                       6567, 13659, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25314, 0, 3,
                                                                       24954, 13389, 24969, 6567,
                                                                       6585, 13689, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25359, 0, 3,
                                                                       24969, 13399, 24984, 6585,
                                                                       6603, 13719, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25404, 0, 3,
                                                                       24999, 13479, 25044, 6639,
                                                                       6675, 13869, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25494, 0, 3,
                                                                       25044, 13509, 25089, 6675,
                                                                       6711, 13929, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25584, 0, 3,
                                                                       25089, 13539, 25134, 6711,
                                                                       6747, 13989, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25674, 0, 3,
                                                                       25134, 13569, 25179, 6747,
                                                                       6783, 14049, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25764, 0, 3,
                                                                       25179, 13599, 25224, 6783,
                                                                       6819, 14109, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25854, 0, 3,
                                                                       25224, 13629, 25269, 6819,
                                                                       6855, 14169, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25944, 0, 3,
                                                                       25269, 13659, 25314, 6855,
                                                                       6891, 14229, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26034, 0, 3,
                                                                       25314, 13689, 25359, 6891,
                                                                       6927, 14289, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26124, 0, 3,
                                                                       25404, 13869, 25494, 6999,
                                                                       7059, 14549, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26274, 0, 3,
                                                                       25494, 13929, 25584, 7059,
                                                                       7119, 14649, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26424, 0, 3,
                                                                       25584, 13989, 25674, 7119,
                                                                       7179, 14749, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26574, 0, 3,
                                                                       25674, 14049, 25764, 7179,
                                                                       7239, 14849, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26724, 0, 3,
                                                                       25764, 14109, 25854, 7239,
                                                                       7299, 14949, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26874, 0, 3,
                                                                       25854, 14169, 25944, 7299,
                                                                       7359, 15049, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27024, 0, 3,
                                                                       25944, 14229, 26034, 7359,
                                                                       7419, 15149, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27174, 0, 3,
                                                                       26124, 14549, 26274, 7539,
                                                                       7629, 15549, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27399, 0, 3,
                                                                       26274, 14649, 26424, 7629,
                                                                       7719, 15699, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27624, 0, 3,
                                                                       26424, 14749, 26574, 7719,
                                                                       7809, 15849, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27849, 0, 3,
                                                                       26574, 14849, 26724, 7809,
                                                                       7899, 15999, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28074, 0, 3,
                                                                       26724, 14949, 26874, 7899,
                                                                       7989, 16149, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28299, 0, 3,
                                                                       26874, 15049, 27024, 7989,
                                                                       8079, 16299, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 28524, 0, 3,
                                                                       27174, 15549, 27399, 8259,
                                                                       8385, 16869, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 28839, 0, 3,
                                                                       27399, 15699, 27624, 8385,
                                                                       8511, 17079, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29154, 0, 3,
                                                                       27624, 15849, 27849, 8511,
                                                                       8637, 17289, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29469, 0, 3,
                                                                       27849, 15999, 28074, 8637,
                                                                       8763, 17499, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29784, 0, 3,
                                                                       28074, 16149, 28299, 8763,
                                                                       8889, 17709, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 30099, 0, 3,
                                                                       28524, 16869, 28839, 9141,
                                                                       9309, 18479, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 30519, 0, 3,
                                                                       28839, 17079, 29154, 9309,
                                                                       9477, 18759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 30939, 0, 3,
                                                                       29154, 17289, 29469, 9477,
                                                                       9645, 19039, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 31359, 0, 3,
                                                                       29469, 17499, 29784, 9645,
                                                                       9813, 19319, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 31779, 0, 3,
                                                                       30099, 18479, 30519,
                                                                       10149, 10365, 20319,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 32319, 0, 3,
                                                                       30519, 18759, 30939,
                                                                       10365, 10581, 20679,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 32859, 0, 3,
                                                                       30939, 19039, 31359,
                                                                       10581, 10797, 21039,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 33399, 0, 3,
                                                                       31779, 20319, 32319,
                                                                       11229, 11499, 22299,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 34074, 0, 3,
                                                                       32319, 20679, 32859,
                                                                       11499, 11769, 22749,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 34749, 0, 3,
                                                                       33399, 22299, 34074,
                                                                       12309, 12639, 24299,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35574, 3, 13299,
                                                                       13309, 24849, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35595, 3, 13309,
                                                                       13319, 24864, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35616, 3, 13319,
                                                                       13329, 24879, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35637, 3, 13329,
                                                                       13339, 24894, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35658, 3, 13339,
                                                                       13349, 24909, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35679, 3, 13349,
                                                                       13359, 24924, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35700, 3, 13359,
                                                                       13369, 24939, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35721, 3, 13369,
                                                                       13379, 24954, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35742, 3, 13379,
                                                                       13389, 24969, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35763, 3, 13389,
                                                                       13399, 24984, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35784, 0, 3,
                                                                       35574, 24849, 35595,
                                                                       13419, 13449, 24999,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35847, 0, 3,
                                                                       35595, 24864, 35616,
                                                                       13449, 13479, 25044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35910, 0, 3,
                                                                       35616, 24879, 35637,
                                                                       13479, 13509, 25089,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35973, 0, 3,
                                                                       35637, 24894, 35658,
                                                                       13509, 13539, 25134,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36036, 0, 3,
                                                                       35658, 24909, 35679,
                                                                       13539, 13569, 25179,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36099, 0, 3,
                                                                       35679, 24924, 35700,
                                                                       13569, 13599, 25224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36162, 0, 3,
                                                                       35700, 24939, 35721,
                                                                       13599, 13629, 25269,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36225, 0, 3,
                                                                       35721, 24954, 35742,
                                                                       13629, 13659, 25314,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36288, 0, 3,
                                                                       35742, 24969, 35763,
                                                                       13659, 13689, 25359,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36351, 0, 3,
                                                                       35784, 24999, 35847,
                                                                       13749, 13809, 25404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36477, 0, 3,
                                                                       35847, 25044, 35910,
                                                                       13809, 13869, 25494,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36603, 0, 3,
                                                                       35910, 25089, 35973,
                                                                       13869, 13929, 25584,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36729, 0, 3,
                                                                       35973, 25134, 36036,
                                                                       13929, 13989, 25674,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36855, 0, 3,
                                                                       36036, 25179, 36099,
                                                                       13989, 14049, 25764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36981, 0, 3,
                                                                       36099, 25224, 36162,
                                                                       14049, 14109, 25854,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37107, 0, 3,
                                                                       36162, 25269, 36225,
                                                                       14109, 14169, 25944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37233, 0, 3,
                                                                       36225, 25314, 36288,
                                                                       14169, 14229, 26034,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37359, 0, 3,
                                                                       36351, 25404, 36477,
                                                                       14349, 14449, 26124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37569, 0, 3,
                                                                       36477, 25494, 36603,
                                                                       14449, 14549, 26274,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37779, 0, 3,
                                                                       36603, 25584, 36729,
                                                                       14549, 14649, 26424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37989, 0, 3,
                                                                       36729, 25674, 36855,
                                                                       14649, 14749, 26574,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38199, 0, 3,
                                                                       36855, 25764, 36981,
                                                                       14749, 14849, 26724,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38409, 0, 3,
                                                                       36981, 25854, 37107,
                                                                       14849, 14949, 26874,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38619, 0, 3,
                                                                       37107, 25944, 37233,
                                                                       14949, 15049, 27024,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 38829, 0, 3,
                                                                       37359, 26124, 37569,
                                                                       15249, 15399, 27174,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 39144, 0, 3,
                                                                       37569, 26274, 37779,
                                                                       15399, 15549, 27399,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 39459, 0, 3,
                                                                       37779, 26424, 37989,
                                                                       15549, 15699, 27624,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 39774, 0, 3,
                                                                       37989, 26574, 38199,
                                                                       15699, 15849, 27849,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40089, 0, 3,
                                                                       38199, 26724, 38409,
                                                                       15849, 15999, 28074,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40404, 0, 3,
                                                                       38409, 26874, 38619,
                                                                       15999, 16149, 28299,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 40719, 0, 3,
                                                                       38829, 27174, 39144,
                                                                       16449, 16659, 28524,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 41160, 0, 3,
                                                                       39144, 27399, 39459,
                                                                       16659, 16869, 28839,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 41601, 0, 3,
                                                                       39459, 27624, 39774,
                                                                       16869, 17079, 29154,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 42042, 0, 3,
                                                                       39774, 27849, 40089,
                                                                       17079, 17289, 29469,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 42483, 0, 3,
                                                                       40089, 28074, 40404,
                                                                       17289, 17499, 29784,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 42924, 0, 3,
                                                                       40719, 28524, 41160,
                                                                       17919, 18199, 30099,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 43512, 0, 3,
                                                                       41160, 28839, 41601,
                                                                       18199, 18479, 30519,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 44100, 0, 3,
                                                                       41601, 29154, 42042,
                                                                       18479, 18759, 30939,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 44688, 0, 3,
                                                                       42042, 29469, 42483,
                                                                       18759, 19039, 31359,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 45276, 0, 3,
                                                                       42924, 30099, 43512,
                                                                       19599, 19959, 31779,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 46032, 0, 3,
                                                                       43512, 30519, 44100,
                                                                       19959, 20319, 32319,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 46788, 0, 3,
                                                                       44100, 30939, 44688,
                                                                       20319, 20679, 32859,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 47544, 0, 3,
                                                                       45276, 31779, 46032,
                                                                       21399, 21849, 33399,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 48489, 0, 3,
                                                                       46032, 32319, 46788,
                                                                       21849, 22299, 34074,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 49434, 0, 3,
                                                                       47544, 33399, 48489,
                                                                       23199, 23749, 34749,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 50589, 40719, 441, ncols);

                    simdfunc::contract_primitives(buffer, 51261, 42924, 588, ncols);

                    simdfunc::contract_primitives(buffer, 52157, 45276, 756, ncols);

                    simdfunc::contract_primitives(buffer, 53309, 47544, 945, ncols);

                    simdfunc::contract_primitives(buffer, 54749, 49434, 1155, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 51030, 50589, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 51849, 51261, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 52913, 52157, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 54254, 53309, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 55904, 54749, 55, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 56509, 51030, 51849, 11, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 57202, 51849, 52913, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 58126, 52913, 54254, 11, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 59314, 54254, 55904, 11, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 60799, 56509, 57202, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 62185, 57202, 58126, 11, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 64033, 58126, 59314, 11, nmax);

        simdtrf::compute_hrr_fh(buffer, coordinates, 66409, 60799, 62185, 11, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 68719, 62185, 64033, 11, nmax);

        simdtrf::compute_hrr_gh(buffer, coordinates, 71799, 66409, 68719, 11, nmax);

        simdtrf::transform_h_inner(buffer, 75264, 71799, 15, 11, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 75264, 121, nmax);
    }

    for (size_t m = 0; m < 1089; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
