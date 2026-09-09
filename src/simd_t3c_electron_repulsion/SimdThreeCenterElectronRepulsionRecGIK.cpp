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


#include "SimdThreeCenterElectronRepulsionRecGIK.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
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
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferDL.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_gik_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_gik_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 213614, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1755 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 213614, 165869, 10740, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16, 17}, ncols, fj, mu, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 66, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 69, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 72, 0, 3, 7, 8,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 78, 0, 3, 8, 9,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 84, 0, 3, 9, 10,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 90, 0, 3, 10, 11,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 96, 0, 3, 11, 12,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 102, 0, 3, 12, 13,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 108, 0, 3, 13, 14,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 114, 0, 3, 14, 15,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 120, 0, 3, 15, 16,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 126, 0, 3, 16, 17,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 132, 0, 3, 17, 18,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 138, 0, 3, 18, 19,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 144, 0, 3, 19, 20,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 150, 0, 3, 20, 21,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 156, 0, 3, 21, 22,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 162, 0, 3, 24, 27,
                                                                       72, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 172, 0, 3, 27, 30,
                                                                       78, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 182, 0, 3, 30, 33,
                                                                       84, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 33, 36,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 202, 0, 3, 36, 39,
                                                                       96, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 212, 0, 3, 39, 42,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 222, 0, 3, 42, 45,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 232, 0, 3, 45, 48,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 242, 0, 3, 48, 51,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 252, 0, 3, 51, 54,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 262, 0, 3, 54, 57,
                                                                       132, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 272, 0, 3, 57, 60,
                                                                       138, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 282, 0, 3, 60, 63,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 292, 0, 3, 63, 66,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 302, 0, 3, 72, 78,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 317, 0, 3, 78, 84,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 332, 0, 3, 84, 90,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 347, 0, 3, 90, 96,
                                                                       192, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 362, 0, 3, 96,
                                                                       102, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 377, 0, 3, 102,
                                                                       108, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 392, 0, 3, 108,
                                                                       114, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 407, 0, 3, 114,
                                                                       120, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 422, 0, 3, 120,
                                                                       126, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 437, 0, 3, 126,
                                                                       132, 252, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 452, 0, 3, 132,
                                                                       138, 262, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 467, 0, 3, 138,
                                                                       144, 272, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 482, 0, 3, 144,
                                                                       150, 282, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 497, 0, 3, 162,
                                                                       172, 302, 317, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 172,
                                                                       182, 317, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 539, 0, 3, 182,
                                                                       192, 332, 347, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 560, 0, 3, 192,
                                                                       202, 347, 362, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 581, 0, 3, 202,
                                                                       212, 362, 377, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 602, 0, 3, 212,
                                                                       222, 377, 392, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 623, 0, 3, 222,
                                                                       232, 392, 407, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 644, 0, 3, 232,
                                                                       242, 407, 422, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 665, 0, 3, 242,
                                                                       252, 422, 437, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 686, 0, 3, 252,
                                                                       262, 437, 452, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 707, 0, 3, 262,
                                                                       272, 452, 467, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 728, 0, 3, 272,
                                                                       282, 467, 482, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 749, 0, 3, 302,
                                                                       317, 497, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 777, 0, 3, 317,
                                                                       332, 518, 539, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 805, 0, 3, 332,
                                                                       347, 539, 560, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 833, 0, 3, 347,
                                                                       362, 560, 581, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 861, 0, 3, 362,
                                                                       377, 581, 602, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 889, 0, 3, 377,
                                                                       392, 602, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 917, 0, 3, 392,
                                                                       407, 623, 644, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 945, 0, 3, 407,
                                                                       422, 644, 665, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 973, 0, 3, 422,
                                                                       437, 665, 686, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1001, 0, 3, 437,
                                                                       452, 686, 707, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1029, 0, 3, 452,
                                                                       467, 707, 728, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1057, 0, 3, 497,
                                                                       518, 749, 777, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1093, 0, 3, 518,
                                                                       539, 777, 805, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1129, 0, 3, 539,
                                                                       560, 805, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1165, 0, 3, 560,
                                                                       581, 833, 861, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1201, 0, 3, 581,
                                                                       602, 861, 889, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1237, 0, 3, 602,
                                                                       623, 889, 917, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1273, 0, 3, 623,
                                                                       644, 917, 945, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1309, 0, 3, 644,
                                                                       665, 945, 973, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1345, 0, 3, 665,
                                                                       686, 973, 1001, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1381, 0, 3, 686,
                                                                       707, 1001, 1029, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1417, 0, 3, 749,
                                                                       777, 1057, 1093, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1462, 0, 3, 777,
                                                                       805, 1093, 1129, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1507, 0, 3, 805,
                                                                       833, 1129, 1165, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1552, 0, 3, 833,
                                                                       861, 1165, 1201, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1597, 0, 3, 861,
                                                                       889, 1201, 1237, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1642, 0, 3, 889,
                                                                       917, 1237, 1273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1687, 0, 3, 917,
                                                                       945, 1273, 1309, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1732, 0, 3, 945,
                                                                       973, 1309, 1345, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1777, 0, 3, 973,
                                                                       1001, 1345, 1381, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1822, 0, 3, 1057,
                                                                       1093, 1417, 1462, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1877, 0, 3, 1093,
                                                                       1129, 1462, 1507, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1932, 0, 3, 1129,
                                                                       1165, 1507, 1552, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1987, 0, 3, 1165,
                                                                       1201, 1552, 1597, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2042, 0, 3, 1201,
                                                                       1237, 1597, 1642, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2097, 0, 3, 1237,
                                                                       1273, 1642, 1687, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2152, 0, 3, 1273,
                                                                       1309, 1687, 1732, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2207, 0, 3, 1309,
                                                                       1345, 1732, 1777, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2262, 0, 3, 1417,
                                                                       1462, 1822, 1877, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2328, 0, 3, 1462,
                                                                       1507, 1877, 1932, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2394, 0, 3, 1507,
                                                                       1552, 1932, 1987, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2460, 0, 3, 1552,
                                                                       1597, 1987, 2042, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2526, 0, 3, 1597,
                                                                       1642, 2042, 2097, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2592, 0, 3, 1642,
                                                                       1687, 2097, 2152, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2658, 0, 3, 1687,
                                                                       1732, 2152, 2207, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2724, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2727, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2730, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2733, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2736, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2739, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2742, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2745, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2748, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2751, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2754, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2757, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2760, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2763, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2766, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2769, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2772, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2775, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2784, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2793, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2802, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2811, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2820, 3, 14, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2829, 3, 15, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2838, 3, 16, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2847, 3, 17, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2856, 3, 18, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2865, 3, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2874, 3, 20, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2883, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2892, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2901, 3, 24, 72,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2919, 3, 27, 78,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2937, 3, 30, 84,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2955, 3, 33, 90,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2973, 3, 36, 96,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2991, 3, 39, 102,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3009, 3, 42, 108,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3027, 3, 45, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3045, 3, 48, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3063, 3, 51, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3081, 3, 54, 132,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3099, 3, 57, 138,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3117, 3, 60, 144,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3135, 3, 63, 150,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3153, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3171, 3, 72, 162,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3201, 3, 78, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3231, 3, 84, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3261, 3, 90, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3291, 3, 96, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3321, 3, 102, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3351, 3, 108, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3381, 3, 114, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3411, 3, 120, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3441, 3, 126, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3471, 3, 132, 262,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3501, 3, 138, 272,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3531, 3, 144, 282,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3561, 3, 150, 292,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3591, 3, 162, 302,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3636, 3, 172, 317,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3681, 3, 182, 332,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3726, 3, 192, 347,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3771, 3, 202, 362,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3816, 3, 212, 377,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3861, 3, 222, 392,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3906, 3, 232, 407,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3951, 3, 242, 422,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3996, 3, 252, 437,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4041, 3, 262, 452,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4086, 3, 272, 467,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4131, 3, 282, 482,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4176, 3, 302, 497,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4239, 3, 317, 518,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4302, 3, 332, 539,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4365, 3, 347, 560,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4428, 3, 362, 581,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4491, 3, 377, 602,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4554, 3, 392, 623,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4617, 3, 407, 644,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4680, 3, 422, 665,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4743, 3, 437, 686,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4806, 3, 452, 707,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4869, 3, 467, 728,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4932, 3, 497, 749,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5016, 3, 518, 777,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5100, 3, 539, 805,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5184, 3, 560, 833,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5268, 3, 581, 861,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5352, 3, 602, 889,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5436, 3, 623, 917,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5520, 3, 644, 945,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5604, 3, 665, 973,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5688, 3, 686,
                                                                       1001, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5772, 3, 707,
                                                                       1029, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5856, 3, 749,
                                                                       1057, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5964, 3, 777,
                                                                       1093, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6072, 3, 805,
                                                                       1129, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6180, 3, 833,
                                                                       1165, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6288, 3, 861,
                                                                       1201, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6396, 3, 889,
                                                                       1237, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6504, 3, 917,
                                                                       1273, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6612, 3, 945,
                                                                       1309, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6720, 3, 973,
                                                                       1345, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6828, 3, 1001,
                                                                       1381, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6936, 3, 1057,
                                                                       1417, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7071, 3, 1093,
                                                                       1462, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7206, 3, 1129,
                                                                       1507, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7341, 3, 1165,
                                                                       1552, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7476, 3, 1201,
                                                                       1597, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7611, 3, 1237,
                                                                       1642, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7746, 3, 1273,
                                                                       1687, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7881, 3, 1309,
                                                                       1732, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8016, 3, 1345,
                                                                       1777, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8151, 3, 1417,
                                                                       1822, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8316, 3, 1462,
                                                                       1877, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8481, 3, 1507,
                                                                       1932, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8646, 3, 1552,
                                                                       1987, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8811, 3, 1597,
                                                                       2042, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8976, 3, 1642,
                                                                       2097, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9141, 3, 1687,
                                                                       2152, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9306, 3, 1732,
                                                                       2207, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9471, 3, 1822,
                                                                       2262, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9669, 3, 1877,
                                                                       2328, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9867, 3, 1932,
                                                                       2394, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10065, 3, 1987,
                                                                       2460, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10263, 3, 2042,
                                                                       2526, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10461, 3, 2097,
                                                                       2592, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10659, 3, 2152,
                                                                       2658, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10857, 3, 7, 8,
                                                                       2730, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10863, 3, 8, 9,
                                                                       2733, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10869, 3, 9, 10,
                                                                       2736, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10875, 3, 10, 11,
                                                                       2739, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10881, 3, 11, 12,
                                                                       2742, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10887, 3, 12, 13,
                                                                       2745, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10893, 3, 13, 14,
                                                                       2748, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10899, 3, 14, 15,
                                                                       2751, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10905, 3, 15, 16,
                                                                       2754, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10911, 3, 16, 17,
                                                                       2757, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10917, 3, 17, 18,
                                                                       2760, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10923, 3, 18, 19,
                                                                       2763, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10929, 3, 19, 20,
                                                                       2766, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10935, 3, 20, 21,
                                                                       2769, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10941, 3, 21, 22,
                                                                       2772, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10947, 0, 3,
                                                                       10857, 2730, 10863, 2775,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10965, 0, 3,
                                                                       10863, 2733, 10869, 2784,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10983, 0, 3,
                                                                       10869, 2736, 10875, 2793,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11001, 0, 3,
                                                                       10875, 2739, 10881, 2802,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11019, 0, 3,
                                                                       10881, 2742, 10887, 2811,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11037, 0, 3,
                                                                       10887, 2745, 10893, 2820,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11055, 0, 3,
                                                                       10893, 2748, 10899, 2829,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11073, 0, 3,
                                                                       10899, 2751, 10905, 2838,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11091, 0, 3,
                                                                       10905, 2754, 10911, 2847,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11109, 0, 3,
                                                                       10911, 2757, 10917, 2856,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11127, 0, 3,
                                                                       10917, 2760, 10923, 2865,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11145, 0, 3,
                                                                       10923, 2763, 10929, 2874,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11163, 0, 3,
                                                                       10929, 2766, 10935, 2883,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11181, 0, 3,
                                                                       10935, 2769, 10941, 2892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11199, 0, 3,
                                                                       10947, 2775, 10965, 72,
                                                                       78, 2937, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11235, 0, 3,
                                                                       10965, 2784, 10983, 78,
                                                                       84, 2955, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11271, 0, 3,
                                                                       10983, 2793, 11001, 84,
                                                                       90, 2973, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11307, 0, 3,
                                                                       11001, 2802, 11019, 90,
                                                                       96, 2991, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11343, 0, 3,
                                                                       11019, 2811, 11037, 96,
                                                                       102, 3009, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11379, 0, 3,
                                                                       11037, 2820, 11055, 102,
                                                                       108, 3027, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11415, 0, 3,
                                                                       11055, 2829, 11073, 108,
                                                                       114, 3045, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11451, 0, 3,
                                                                       11073, 2838, 11091, 114,
                                                                       120, 3063, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11487, 0, 3,
                                                                       11091, 2847, 11109, 120,
                                                                       126, 3081, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11523, 0, 3,
                                                                       11109, 2856, 11127, 126,
                                                                       132, 3099, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11559, 0, 3,
                                                                       11127, 2865, 11145, 132,
                                                                       138, 3117, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11595, 0, 3,
                                                                       11145, 2874, 11163, 138,
                                                                       144, 3135, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11631, 0, 3,
                                                                       11163, 2883, 11181, 144,
                                                                       150, 3153, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11667, 0, 3,
                                                                       11199, 2937, 11235, 162,
                                                                       172, 3231, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11727, 0, 3,
                                                                       11235, 2955, 11271, 172,
                                                                       182, 3261, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11787, 0, 3,
                                                                       11271, 2973, 11307, 182,
                                                                       192, 3291, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11847, 0, 3,
                                                                       11307, 2991, 11343, 192,
                                                                       202, 3321, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11907, 0, 3,
                                                                       11343, 3009, 11379, 202,
                                                                       212, 3351, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11967, 0, 3,
                                                                       11379, 3027, 11415, 212,
                                                                       222, 3381, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12027, 0, 3,
                                                                       11415, 3045, 11451, 222,
                                                                       232, 3411, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12087, 0, 3,
                                                                       11451, 3063, 11487, 232,
                                                                       242, 3441, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12147, 0, 3,
                                                                       11487, 3081, 11523, 242,
                                                                       252, 3471, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12207, 0, 3,
                                                                       11523, 3099, 11559, 252,
                                                                       262, 3501, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12267, 0, 3,
                                                                       11559, 3117, 11595, 262,
                                                                       272, 3531, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12327, 0, 3,
                                                                       11595, 3135, 11631, 272,
                                                                       282, 3561, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12387, 0, 3,
                                                                       11667, 3231, 11727, 302,
                                                                       317, 3681, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12477, 0, 3,
                                                                       11727, 3261, 11787, 317,
                                                                       332, 3726, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12567, 0, 3,
                                                                       11787, 3291, 11847, 332,
                                                                       347, 3771, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12657, 0, 3,
                                                                       11847, 3321, 11907, 347,
                                                                       362, 3816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12747, 0, 3,
                                                                       11907, 3351, 11967, 362,
                                                                       377, 3861, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12837, 0, 3,
                                                                       11967, 3381, 12027, 377,
                                                                       392, 3906, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12927, 0, 3,
                                                                       12027, 3411, 12087, 392,
                                                                       407, 3951, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13017, 0, 3,
                                                                       12087, 3441, 12147, 407,
                                                                       422, 3996, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13107, 0, 3,
                                                                       12147, 3471, 12207, 422,
                                                                       437, 4041, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13197, 0, 3,
                                                                       12207, 3501, 12267, 437,
                                                                       452, 4086, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13287, 0, 3,
                                                                       12267, 3531, 12327, 452,
                                                                       467, 4131, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13377, 0, 3,
                                                                       12387, 3681, 12477, 497,
                                                                       518, 4302, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13503, 0, 3,
                                                                       12477, 3726, 12567, 518,
                                                                       539, 4365, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13629, 0, 3,
                                                                       12567, 3771, 12657, 539,
                                                                       560, 4428, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13755, 0, 3,
                                                                       12657, 3816, 12747, 560,
                                                                       581, 4491, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13881, 0, 3,
                                                                       12747, 3861, 12837, 581,
                                                                       602, 4554, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14007, 0, 3,
                                                                       12837, 3906, 12927, 602,
                                                                       623, 4617, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14133, 0, 3,
                                                                       12927, 3951, 13017, 623,
                                                                       644, 4680, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14259, 0, 3,
                                                                       13017, 3996, 13107, 644,
                                                                       665, 4743, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14385, 0, 3,
                                                                       13107, 4041, 13197, 665,
                                                                       686, 4806, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14511, 0, 3,
                                                                       13197, 4086, 13287, 686,
                                                                       707, 4869, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14637, 0, 3,
                                                                       13377, 4302, 13503, 749,
                                                                       777, 5100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14805, 0, 3,
                                                                       13503, 4365, 13629, 777,
                                                                       805, 5184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14973, 0, 3,
                                                                       13629, 4428, 13755, 805,
                                                                       833, 5268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15141, 0, 3,
                                                                       13755, 4491, 13881, 833,
                                                                       861, 5352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15309, 0, 3,
                                                                       13881, 4554, 14007, 861,
                                                                       889, 5436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15477, 0, 3,
                                                                       14007, 4617, 14133, 889,
                                                                       917, 5520, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15645, 0, 3,
                                                                       14133, 4680, 14259, 917,
                                                                       945, 5604, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15813, 0, 3,
                                                                       14259, 4743, 14385, 945,
                                                                       973, 5688, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15981, 0, 3,
                                                                       14385, 4806, 14511, 973,
                                                                       1001, 5772, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16149, 0, 3,
                                                                       14637, 5100, 14805, 1057,
                                                                       1093, 6072, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16365, 0, 3,
                                                                       14805, 5184, 14973, 1093,
                                                                       1129, 6180, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16581, 0, 3,
                                                                       14973, 5268, 15141, 1129,
                                                                       1165, 6288, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16797, 0, 3,
                                                                       15141, 5352, 15309, 1165,
                                                                       1201, 6396, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17013, 0, 3,
                                                                       15309, 5436, 15477, 1201,
                                                                       1237, 6504, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17229, 0, 3,
                                                                       15477, 5520, 15645, 1237,
                                                                       1273, 6612, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17445, 0, 3,
                                                                       15645, 5604, 15813, 1273,
                                                                       1309, 6720, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17661, 0, 3,
                                                                       15813, 5688, 15981, 1309,
                                                                       1345, 6828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17877, 0, 3,
                                                                       16149, 6072, 16365, 1417,
                                                                       1462, 7206, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18147, 0, 3,
                                                                       16365, 6180, 16581, 1462,
                                                                       1507, 7341, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18417, 0, 3,
                                                                       16581, 6288, 16797, 1507,
                                                                       1552, 7476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18687, 0, 3,
                                                                       16797, 6396, 17013, 1552,
                                                                       1597, 7611, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18957, 0, 3,
                                                                       17013, 6504, 17229, 1597,
                                                                       1642, 7746, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 19227, 0, 3,
                                                                       17229, 6612, 17445, 1642,
                                                                       1687, 7881, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 19497, 0, 3,
                                                                       17445, 6720, 17661, 1687,
                                                                       1732, 8016, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19767, 0, 3,
                                                                       17877, 7206, 18147, 1822,
                                                                       1877, 8481, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20097, 0, 3,
                                                                       18147, 7341, 18417, 1877,
                                                                       1932, 8646, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20427, 0, 3,
                                                                       18417, 7476, 18687, 1932,
                                                                       1987, 8811, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20757, 0, 3,
                                                                       18687, 7611, 18957, 1987,
                                                                       2042, 8976, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 21087, 0, 3,
                                                                       18957, 7746, 19227, 2042,
                                                                       2097, 9141, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 21417, 0, 3,
                                                                       19227, 7881, 19497, 2097,
                                                                       2152, 9306, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 21747, 0, 3,
                                                                       19767, 8481, 20097, 2262,
                                                                       2328, 9867, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 22143, 0, 3,
                                                                       20097, 8646, 20427, 2328,
                                                                       2394, 10065, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 22539, 0, 3,
                                                                       20427, 8811, 20757, 2394,
                                                                       2460, 10263, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 22935, 0, 3,
                                                                       20757, 8976, 21087, 2460,
                                                                       2526, 10461, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 23331, 0, 3,
                                                                       21087, 9141, 21417, 2526,
                                                                       2592, 10659, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23727, 3, 2724,
                                                                       2727, 10857, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23737, 3, 2727,
                                                                       2730, 10863, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23747, 3, 2730,
                                                                       2733, 10869, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23757, 3, 2733,
                                                                       2736, 10875, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23767, 3, 2736,
                                                                       2739, 10881, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23777, 3, 2739,
                                                                       2742, 10887, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23787, 3, 2742,
                                                                       2745, 10893, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23797, 3, 2745,
                                                                       2748, 10899, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23807, 3, 2748,
                                                                       2751, 10905, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23817, 3, 2751,
                                                                       2754, 10911, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23827, 3, 2754,
                                                                       2757, 10917, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23837, 3, 2757,
                                                                       2760, 10923, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23847, 3, 2760,
                                                                       2763, 10929, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23857, 3, 2763,
                                                                       2766, 10935, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23867, 3, 2766,
                                                                       2769, 10941, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23877, 0, 3,
                                                                       23727, 10857, 23737,
                                                                       10947, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23907, 0, 3,
                                                                       23737, 10863, 23747,
                                                                       10965, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23937, 0, 3,
                                                                       23747, 10869, 23757,
                                                                       10983, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23967, 0, 3,
                                                                       23757, 10875, 23767,
                                                                       11001, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23997, 0, 3,
                                                                       23767, 10881, 23777,
                                                                       11019, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24027, 0, 3,
                                                                       23777, 10887, 23787,
                                                                       11037, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24057, 0, 3,
                                                                       23787, 10893, 23797,
                                                                       11055, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24087, 0, 3,
                                                                       23797, 10899, 23807,
                                                                       11073, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24117, 0, 3,
                                                                       23807, 10905, 23817,
                                                                       11091, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24147, 0, 3,
                                                                       23817, 10911, 23827,
                                                                       11109, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24177, 0, 3,
                                                                       23827, 10917, 23837,
                                                                       11127, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24207, 0, 3,
                                                                       23837, 10923, 23847,
                                                                       11145, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24237, 0, 3,
                                                                       23847, 10929, 23857,
                                                                       11163, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24267, 0, 3,
                                                                       23857, 10935, 23867,
                                                                       11181, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24297, 0, 3,
                                                                       23877, 10947, 23907, 2901,
                                                                       2919, 11199, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24357, 0, 3,
                                                                       23907, 10965, 23937, 2919,
                                                                       2937, 11235, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24417, 0, 3,
                                                                       23937, 10983, 23967, 2937,
                                                                       2955, 11271, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24477, 0, 3,
                                                                       23967, 11001, 23997, 2955,
                                                                       2973, 11307, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24537, 0, 3,
                                                                       23997, 11019, 24027, 2973,
                                                                       2991, 11343, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24597, 0, 3,
                                                                       24027, 11037, 24057, 2991,
                                                                       3009, 11379, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24657, 0, 3,
                                                                       24057, 11055, 24087, 3009,
                                                                       3027, 11415, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24717, 0, 3,
                                                                       24087, 11073, 24117, 3027,
                                                                       3045, 11451, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24777, 0, 3,
                                                                       24117, 11091, 24147, 3045,
                                                                       3063, 11487, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24837, 0, 3,
                                                                       24147, 11109, 24177, 3063,
                                                                       3081, 11523, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24897, 0, 3,
                                                                       24177, 11127, 24207, 3081,
                                                                       3099, 11559, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24957, 0, 3,
                                                                       24207, 11145, 24237, 3099,
                                                                       3117, 11595, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25017, 0, 3,
                                                                       24237, 11163, 24267, 3117,
                                                                       3135, 11631, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25077, 0, 3,
                                                                       24297, 11199, 24357, 3171,
                                                                       3201, 11667, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25177, 0, 3,
                                                                       24357, 11235, 24417, 3201,
                                                                       3231, 11727, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25277, 0, 3,
                                                                       24417, 11271, 24477, 3231,
                                                                       3261, 11787, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25377, 0, 3,
                                                                       24477, 11307, 24537, 3261,
                                                                       3291, 11847, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25477, 0, 3,
                                                                       24537, 11343, 24597, 3291,
                                                                       3321, 11907, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25577, 0, 3,
                                                                       24597, 11379, 24657, 3321,
                                                                       3351, 11967, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25677, 0, 3,
                                                                       24657, 11415, 24717, 3351,
                                                                       3381, 12027, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25777, 0, 3,
                                                                       24717, 11451, 24777, 3381,
                                                                       3411, 12087, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25877, 0, 3,
                                                                       24777, 11487, 24837, 3411,
                                                                       3441, 12147, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25977, 0, 3,
                                                                       24837, 11523, 24897, 3441,
                                                                       3471, 12207, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26077, 0, 3,
                                                                       24897, 11559, 24957, 3471,
                                                                       3501, 12267, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26177, 0, 3,
                                                                       24957, 11595, 25017, 3501,
                                                                       3531, 12327, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26277, 0, 3,
                                                                       25077, 11667, 25177, 3591,
                                                                       3636, 12387, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26427, 0, 3,
                                                                       25177, 11727, 25277, 3636,
                                                                       3681, 12477, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26577, 0, 3,
                                                                       25277, 11787, 25377, 3681,
                                                                       3726, 12567, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26727, 0, 3,
                                                                       25377, 11847, 25477, 3726,
                                                                       3771, 12657, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26877, 0, 3,
                                                                       25477, 11907, 25577, 3771,
                                                                       3816, 12747, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27027, 0, 3,
                                                                       25577, 11967, 25677, 3816,
                                                                       3861, 12837, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27177, 0, 3,
                                                                       25677, 12027, 25777, 3861,
                                                                       3906, 12927, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27327, 0, 3,
                                                                       25777, 12087, 25877, 3906,
                                                                       3951, 13017, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27477, 0, 3,
                                                                       25877, 12147, 25977, 3951,
                                                                       3996, 13107, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27627, 0, 3,
                                                                       25977, 12207, 26077, 3996,
                                                                       4041, 13197, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27777, 0, 3,
                                                                       26077, 12267, 26177, 4041,
                                                                       4086, 13287, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 27927, 0, 3,
                                                                       26277, 12387, 26427, 4176,
                                                                       4239, 13377, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28137, 0, 3,
                                                                       26427, 12477, 26577, 4239,
                                                                       4302, 13503, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28347, 0, 3,
                                                                       26577, 12567, 26727, 4302,
                                                                       4365, 13629, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28557, 0, 3,
                                                                       26727, 12657, 26877, 4365,
                                                                       4428, 13755, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28767, 0, 3,
                                                                       26877, 12747, 27027, 4428,
                                                                       4491, 13881, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28977, 0, 3,
                                                                       27027, 12837, 27177, 4491,
                                                                       4554, 14007, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29187, 0, 3,
                                                                       27177, 12927, 27327, 4554,
                                                                       4617, 14133, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29397, 0, 3,
                                                                       27327, 13017, 27477, 4617,
                                                                       4680, 14259, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29607, 0, 3,
                                                                       27477, 13107, 27627, 4680,
                                                                       4743, 14385, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29817, 0, 3,
                                                                       27627, 13197, 27777, 4743,
                                                                       4806, 14511, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30027, 0, 3,
                                                                       27927, 13377, 28137, 4932,
                                                                       5016, 14637, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30307, 0, 3,
                                                                       28137, 13503, 28347, 5016,
                                                                       5100, 14805, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30587, 0, 3,
                                                                       28347, 13629, 28557, 5100,
                                                                       5184, 14973, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30867, 0, 3,
                                                                       28557, 13755, 28767, 5184,
                                                                       5268, 15141, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31147, 0, 3,
                                                                       28767, 13881, 28977, 5268,
                                                                       5352, 15309, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31427, 0, 3,
                                                                       28977, 14007, 29187, 5352,
                                                                       5436, 15477, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31707, 0, 3,
                                                                       29187, 14133, 29397, 5436,
                                                                       5520, 15645, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31987, 0, 3,
                                                                       29397, 14259, 29607, 5520,
                                                                       5604, 15813, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32267, 0, 3,
                                                                       29607, 14385, 29817, 5604,
                                                                       5688, 15981, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32547, 0, 3,
                                                                       30027, 14637, 30307, 5856,
                                                                       5964, 16149, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32907, 0, 3,
                                                                       30307, 14805, 30587, 5964,
                                                                       6072, 16365, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33267, 0, 3,
                                                                       30587, 14973, 30867, 6072,
                                                                       6180, 16581, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33627, 0, 3,
                                                                       30867, 15141, 31147, 6180,
                                                                       6288, 16797, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33987, 0, 3,
                                                                       31147, 15309, 31427, 6288,
                                                                       6396, 17013, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 34347, 0, 3,
                                                                       31427, 15477, 31707, 6396,
                                                                       6504, 17229, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 34707, 0, 3,
                                                                       31707, 15645, 31987, 6504,
                                                                       6612, 17445, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 35067, 0, 3,
                                                                       31987, 15813, 32267, 6612,
                                                                       6720, 17661, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 35427, 0, 3,
                                                                       32547, 16149, 32907, 6936,
                                                                       7071, 17877, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 35877, 0, 3,
                                                                       32907, 16365, 33267, 7071,
                                                                       7206, 18147, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 36327, 0, 3,
                                                                       33267, 16581, 33627, 7206,
                                                                       7341, 18417, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 36777, 0, 3,
                                                                       33627, 16797, 33987, 7341,
                                                                       7476, 18687, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 37227, 0, 3,
                                                                       33987, 17013, 34347, 7476,
                                                                       7611, 18957, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 37677, 0, 3,
                                                                       34347, 17229, 34707, 7611,
                                                                       7746, 19227, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 38127, 0, 3,
                                                                       34707, 17445, 35067, 7746,
                                                                       7881, 19497, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 38577, 0, 3,
                                                                       35427, 17877, 35877, 8151,
                                                                       8316, 19767, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 39127, 0, 3,
                                                                       35877, 18147, 36327, 8316,
                                                                       8481, 20097, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 39677, 0, 3,
                                                                       36327, 18417, 36777, 8481,
                                                                       8646, 20427, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 40227, 0, 3,
                                                                       36777, 18687, 37227, 8646,
                                                                       8811, 20757, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 40777, 0, 3,
                                                                       37227, 18957, 37677, 8811,
                                                                       8976, 21087, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 41327, 0, 3,
                                                                       37677, 19227, 38127, 8976,
                                                                       9141, 21417, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 41877, 0, 3,
                                                                       38577, 19767, 39127, 9471,
                                                                       9669, 21747, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 42537, 0, 3,
                                                                       39127, 20097, 39677, 9669,
                                                                       9867, 22143, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 43197, 0, 3,
                                                                       39677, 20427, 40227, 9867,
                                                                       10065, 22539, ncols,
                                                                       gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 43857, 0, 3,
                                                                       40227, 20757, 40777,
                                                                       10065, 10263, 22935,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 44517, 0, 3,
                                                                       40777, 21087, 41327,
                                                                       10263, 10461, 23331,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45177, 3, 10857,
                                                                       10863, 23747, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45192, 3, 10863,
                                                                       10869, 23757, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45207, 3, 10869,
                                                                       10875, 23767, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45222, 3, 10875,
                                                                       10881, 23777, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45237, 3, 10881,
                                                                       10887, 23787, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45252, 3, 10887,
                                                                       10893, 23797, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45267, 3, 10893,
                                                                       10899, 23807, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45282, 3, 10899,
                                                                       10905, 23817, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45297, 3, 10905,
                                                                       10911, 23827, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45312, 3, 10911,
                                                                       10917, 23837, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45327, 3, 10917,
                                                                       10923, 23847, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45342, 3, 10923,
                                                                       10929, 23857, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45357, 3, 10929,
                                                                       10935, 23867, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45372, 0, 3,
                                                                       45177, 23747, 45192,
                                                                       10947, 10965, 23937,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45417, 0, 3,
                                                                       45192, 23757, 45207,
                                                                       10965, 10983, 23967,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45462, 0, 3,
                                                                       45207, 23767, 45222,
                                                                       10983, 11001, 23997,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45507, 0, 3,
                                                                       45222, 23777, 45237,
                                                                       11001, 11019, 24027,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45552, 0, 3,
                                                                       45237, 23787, 45252,
                                                                       11019, 11037, 24057,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45597, 0, 3,
                                                                       45252, 23797, 45267,
                                                                       11037, 11055, 24087,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45642, 0, 3,
                                                                       45267, 23807, 45282,
                                                                       11055, 11073, 24117,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45687, 0, 3,
                                                                       45282, 23817, 45297,
                                                                       11073, 11091, 24147,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45732, 0, 3,
                                                                       45297, 23827, 45312,
                                                                       11091, 11109, 24177,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45777, 0, 3,
                                                                       45312, 23837, 45327,
                                                                       11109, 11127, 24207,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45822, 0, 3,
                                                                       45327, 23847, 45342,
                                                                       11127, 11145, 24237,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45867, 0, 3,
                                                                       45342, 23857, 45357,
                                                                       11145, 11163, 24267,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45912, 0, 3,
                                                                       45372, 23937, 45417,
                                                                       11199, 11235, 24417,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46002, 0, 3,
                                                                       45417, 23967, 45462,
                                                                       11235, 11271, 24477,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46092, 0, 3,
                                                                       45462, 23997, 45507,
                                                                       11271, 11307, 24537,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46182, 0, 3,
                                                                       45507, 24027, 45552,
                                                                       11307, 11343, 24597,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46272, 0, 3,
                                                                       45552, 24057, 45597,
                                                                       11343, 11379, 24657,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46362, 0, 3,
                                                                       45597, 24087, 45642,
                                                                       11379, 11415, 24717,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46452, 0, 3,
                                                                       45642, 24117, 45687,
                                                                       11415, 11451, 24777,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46542, 0, 3,
                                                                       45687, 24147, 45732,
                                                                       11451, 11487, 24837,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46632, 0, 3,
                                                                       45732, 24177, 45777,
                                                                       11487, 11523, 24897,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46722, 0, 3,
                                                                       45777, 24207, 45822,
                                                                       11523, 11559, 24957,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46812, 0, 3,
                                                                       45822, 24237, 45867,
                                                                       11559, 11595, 25017,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46902, 0, 3,
                                                                       45912, 24417, 46002,
                                                                       11667, 11727, 25277,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47052, 0, 3,
                                                                       46002, 24477, 46092,
                                                                       11727, 11787, 25377,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47202, 0, 3,
                                                                       46092, 24537, 46182,
                                                                       11787, 11847, 25477,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47352, 0, 3,
                                                                       46182, 24597, 46272,
                                                                       11847, 11907, 25577,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47502, 0, 3,
                                                                       46272, 24657, 46362,
                                                                       11907, 11967, 25677,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47652, 0, 3,
                                                                       46362, 24717, 46452,
                                                                       11967, 12027, 25777,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47802, 0, 3,
                                                                       46452, 24777, 46542,
                                                                       12027, 12087, 25877,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47952, 0, 3,
                                                                       46542, 24837, 46632,
                                                                       12087, 12147, 25977,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 48102, 0, 3,
                                                                       46632, 24897, 46722,
                                                                       12147, 12207, 26077,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 48252, 0, 3,
                                                                       46722, 24957, 46812,
                                                                       12207, 12267, 26177,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48402, 0, 3,
                                                                       46902, 25277, 47052,
                                                                       12387, 12477, 26577,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48627, 0, 3,
                                                                       47052, 25377, 47202,
                                                                       12477, 12567, 26727,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48852, 0, 3,
                                                                       47202, 25477, 47352,
                                                                       12567, 12657, 26877,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49077, 0, 3,
                                                                       47352, 25577, 47502,
                                                                       12657, 12747, 27027,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49302, 0, 3,
                                                                       47502, 25677, 47652,
                                                                       12747, 12837, 27177,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49527, 0, 3,
                                                                       47652, 25777, 47802,
                                                                       12837, 12927, 27327,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49752, 0, 3,
                                                                       47802, 25877, 47952,
                                                                       12927, 13017, 27477,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49977, 0, 3,
                                                                       47952, 25977, 48102,
                                                                       13017, 13107, 27627,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 50202, 0, 3,
                                                                       48102, 26077, 48252,
                                                                       13107, 13197, 27777,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50427, 0, 3,
                                                                       48402, 26577, 48627,
                                                                       13377, 13503, 28347,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50742, 0, 3,
                                                                       48627, 26727, 48852,
                                                                       13503, 13629, 28557,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 51057, 0, 3,
                                                                       48852, 26877, 49077,
                                                                       13629, 13755, 28767,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 51372, 0, 3,
                                                                       49077, 27027, 49302,
                                                                       13755, 13881, 28977,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 51687, 0, 3,
                                                                       49302, 27177, 49527,
                                                                       13881, 14007, 29187,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 52002, 0, 3,
                                                                       49527, 27327, 49752,
                                                                       14007, 14133, 29397,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 52317, 0, 3,
                                                                       49752, 27477, 49977,
                                                                       14133, 14259, 29607,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 52632, 0, 3,
                                                                       49977, 27627, 50202,
                                                                       14259, 14385, 29817,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 52947, 0, 3,
                                                                       50427, 28347, 50742,
                                                                       14637, 14805, 30587,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 53367, 0, 3,
                                                                       50742, 28557, 51057,
                                                                       14805, 14973, 30867,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 53787, 0, 3,
                                                                       51057, 28767, 51372,
                                                                       14973, 15141, 31147,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 54207, 0, 3,
                                                                       51372, 28977, 51687,
                                                                       15141, 15309, 31427,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 54627, 0, 3,
                                                                       51687, 29187, 52002,
                                                                       15309, 15477, 31707,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 55047, 0, 3,
                                                                       52002, 29397, 52317,
                                                                       15477, 15645, 31987,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 55467, 0, 3,
                                                                       52317, 29607, 52632,
                                                                       15645, 15813, 32267,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 55887, 0, 3,
                                                                       52947, 30587, 53367,
                                                                       16149, 16365, 33267,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 56427, 0, 3,
                                                                       53367, 30867, 53787,
                                                                       16365, 16581, 33627,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 56967, 0, 3,
                                                                       53787, 31147, 54207,
                                                                       16581, 16797, 33987,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 57507, 0, 3,
                                                                       54207, 31427, 54627,
                                                                       16797, 17013, 34347,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 58047, 0, 3,
                                                                       54627, 31707, 55047,
                                                                       17013, 17229, 34707,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 58587, 0, 3,
                                                                       55047, 31987, 55467,
                                                                       17229, 17445, 35067,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 59127, 0, 3,
                                                                       55887, 33267, 56427,
                                                                       17877, 18147, 36327,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 59802, 0, 3,
                                                                       56427, 33627, 56967,
                                                                       18147, 18417, 36777,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 60477, 0, 3,
                                                                       56967, 33987, 57507,
                                                                       18417, 18687, 37227,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 61152, 0, 3,
                                                                       57507, 34347, 58047,
                                                                       18687, 18957, 37677,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 61827, 0, 3,
                                                                       58047, 34707, 58587,
                                                                       18957, 19227, 38127,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 62502, 0, 3,
                                                                       59127, 36327, 59802,
                                                                       19767, 20097, 39677,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 63327, 0, 3,
                                                                       59802, 36777, 60477,
                                                                       20097, 20427, 40227,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 64152, 0, 3,
                                                                       60477, 37227, 61152,
                                                                       20427, 20757, 40777,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 64977, 0, 3,
                                                                       61152, 37677, 61827,
                                                                       20757, 21087, 41327,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 65802, 0, 3,
                                                                       62502, 39677, 63327,
                                                                       21747, 22143, 43197,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 66792, 0, 3,
                                                                       63327, 40227, 64152,
                                                                       22143, 22539, 43857,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 67782, 0, 3,
                                                                       64152, 40777, 64977,
                                                                       22539, 22935, 44517,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68772, 3, 23727,
                                                                       23737, 45177, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68793, 3, 23737,
                                                                       23747, 45192, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68814, 3, 23747,
                                                                       23757, 45207, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68835, 3, 23757,
                                                                       23767, 45222, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68856, 3, 23767,
                                                                       23777, 45237, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68877, 3, 23777,
                                                                       23787, 45252, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68898, 3, 23787,
                                                                       23797, 45267, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68919, 3, 23797,
                                                                       23807, 45282, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68940, 3, 23807,
                                                                       23817, 45297, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68961, 3, 23817,
                                                                       23827, 45312, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68982, 3, 23827,
                                                                       23837, 45327, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 69003, 3, 23837,
                                                                       23847, 45342, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 69024, 3, 23847,
                                                                       23857, 45357, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69045, 0, 3,
                                                                       68772, 45177, 68793,
                                                                       23877, 23907, 45372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69108, 0, 3,
                                                                       68793, 45192, 68814,
                                                                       23907, 23937, 45417,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69171, 0, 3,
                                                                       68814, 45207, 68835,
                                                                       23937, 23967, 45462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69234, 0, 3,
                                                                       68835, 45222, 68856,
                                                                       23967, 23997, 45507,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69297, 0, 3,
                                                                       68856, 45237, 68877,
                                                                       23997, 24027, 45552,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69360, 0, 3,
                                                                       68877, 45252, 68898,
                                                                       24027, 24057, 45597,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69423, 0, 3,
                                                                       68898, 45267, 68919,
                                                                       24057, 24087, 45642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69486, 0, 3,
                                                                       68919, 45282, 68940,
                                                                       24087, 24117, 45687,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69549, 0, 3,
                                                                       68940, 45297, 68961,
                                                                       24117, 24147, 45732,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69612, 0, 3,
                                                                       68961, 45312, 68982,
                                                                       24147, 24177, 45777,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69675, 0, 3,
                                                                       68982, 45327, 69003,
                                                                       24177, 24207, 45822,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69738, 0, 3,
                                                                       69003, 45342, 69024,
                                                                       24207, 24237, 45867,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 69801, 0, 3,
                                                                       69045, 45372, 69108,
                                                                       24297, 24357, 45912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 69927, 0, 3,
                                                                       69108, 45417, 69171,
                                                                       24357, 24417, 46002,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70053, 0, 3,
                                                                       69171, 45462, 69234,
                                                                       24417, 24477, 46092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70179, 0, 3,
                                                                       69234, 45507, 69297,
                                                                       24477, 24537, 46182,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70305, 0, 3,
                                                                       69297, 45552, 69360,
                                                                       24537, 24597, 46272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70431, 0, 3,
                                                                       69360, 45597, 69423,
                                                                       24597, 24657, 46362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70557, 0, 3,
                                                                       69423, 45642, 69486,
                                                                       24657, 24717, 46452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70683, 0, 3,
                                                                       69486, 45687, 69549,
                                                                       24717, 24777, 46542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70809, 0, 3,
                                                                       69549, 45732, 69612,
                                                                       24777, 24837, 46632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70935, 0, 3,
                                                                       69612, 45777, 69675,
                                                                       24837, 24897, 46722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 71061, 0, 3,
                                                                       69675, 45822, 69738,
                                                                       24897, 24957, 46812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 71187, 0, 3,
                                                                       69801, 45912, 69927,
                                                                       25077, 25177, 46902,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 71397, 0, 3,
                                                                       69927, 46002, 70053,
                                                                       25177, 25277, 47052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 71607, 0, 3,
                                                                       70053, 46092, 70179,
                                                                       25277, 25377, 47202,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 71817, 0, 3,
                                                                       70179, 46182, 70305,
                                                                       25377, 25477, 47352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 72027, 0, 3,
                                                                       70305, 46272, 70431,
                                                                       25477, 25577, 47502,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 72237, 0, 3,
                                                                       70431, 46362, 70557,
                                                                       25577, 25677, 47652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 72447, 0, 3,
                                                                       70557, 46452, 70683,
                                                                       25677, 25777, 47802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 72657, 0, 3,
                                                                       70683, 46542, 70809,
                                                                       25777, 25877, 47952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 72867, 0, 3,
                                                                       70809, 46632, 70935,
                                                                       25877, 25977, 48102,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 73077, 0, 3,
                                                                       70935, 46722, 71061,
                                                                       25977, 26077, 48252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 73287, 0, 3,
                                                                       71187, 46902, 71397,
                                                                       26277, 26427, 48402,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 73602, 0, 3,
                                                                       71397, 47052, 71607,
                                                                       26427, 26577, 48627,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 73917, 0, 3,
                                                                       71607, 47202, 71817,
                                                                       26577, 26727, 48852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 74232, 0, 3,
                                                                       71817, 47352, 72027,
                                                                       26727, 26877, 49077,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 74547, 0, 3,
                                                                       72027, 47502, 72237,
                                                                       26877, 27027, 49302,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 74862, 0, 3,
                                                                       72237, 47652, 72447,
                                                                       27027, 27177, 49527,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 75177, 0, 3,
                                                                       72447, 47802, 72657,
                                                                       27177, 27327, 49752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 75492, 0, 3,
                                                                       72657, 47952, 72867,
                                                                       27327, 27477, 49977,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 75807, 0, 3,
                                                                       72867, 48102, 73077,
                                                                       27477, 27627, 50202,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 76122, 0, 3,
                                                                       73287, 48402, 73602,
                                                                       27927, 28137, 50427,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 76563, 0, 3,
                                                                       73602, 48627, 73917,
                                                                       28137, 28347, 50742,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 77004, 0, 3,
                                                                       73917, 48852, 74232,
                                                                       28347, 28557, 51057,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 77445, 0, 3,
                                                                       74232, 49077, 74547,
                                                                       28557, 28767, 51372,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 77886, 0, 3,
                                                                       74547, 49302, 74862,
                                                                       28767, 28977, 51687,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 78327, 0, 3,
                                                                       74862, 49527, 75177,
                                                                       28977, 29187, 52002,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 78768, 0, 3,
                                                                       75177, 49752, 75492,
                                                                       29187, 29397, 52317,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 79209, 0, 3,
                                                                       75492, 49977, 75807,
                                                                       29397, 29607, 52632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 79650, 0, 3,
                                                                       76122, 50427, 76563,
                                                                       30027, 30307, 52947,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 80238, 0, 3,
                                                                       76563, 50742, 77004,
                                                                       30307, 30587, 53367,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 80826, 0, 3,
                                                                       77004, 51057, 77445,
                                                                       30587, 30867, 53787,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 81414, 0, 3,
                                                                       77445, 51372, 77886,
                                                                       30867, 31147, 54207,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 82002, 0, 3,
                                                                       77886, 51687, 78327,
                                                                       31147, 31427, 54627,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 82590, 0, 3,
                                                                       78327, 52002, 78768,
                                                                       31427, 31707, 55047,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 83178, 0, 3,
                                                                       78768, 52317, 79209,
                                                                       31707, 31987, 55467,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 83766, 0, 3,
                                                                       79650, 52947, 80238,
                                                                       32547, 32907, 55887,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 84522, 0, 3,
                                                                       80238, 53367, 80826,
                                                                       32907, 33267, 56427,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 85278, 0, 3,
                                                                       80826, 53787, 81414,
                                                                       33267, 33627, 56967,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 86034, 0, 3,
                                                                       81414, 54207, 82002,
                                                                       33627, 33987, 57507,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 86790, 0, 3,
                                                                       82002, 54627, 82590,
                                                                       33987, 34347, 58047,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 87546, 0, 3,
                                                                       82590, 55047, 83178,
                                                                       34347, 34707, 58587,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 88302, 0, 3,
                                                                       83766, 55887, 84522,
                                                                       35427, 35877, 59127,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 89247, 0, 3,
                                                                       84522, 56427, 85278,
                                                                       35877, 36327, 59802,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 90192, 0, 3,
                                                                       85278, 56967, 86034,
                                                                       36327, 36777, 60477,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 91137, 0, 3,
                                                                       86034, 57507, 86790,
                                                                       36777, 37227, 61152,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 92082, 0, 3,
                                                                       86790, 58047, 87546,
                                                                       37227, 37677, 61827,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 93027, 0, 3,
                                                                       88302, 59127, 89247,
                                                                       38577, 39127, 62502,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 94182, 0, 3,
                                                                       89247, 59802, 90192,
                                                                       39127, 39677, 63327,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 95337, 0, 3,
                                                                       90192, 60477, 91137,
                                                                       39677, 40227, 64152,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 96492, 0, 3,
                                                                       91137, 61152, 92082,
                                                                       40227, 40777, 64977,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 97647, 0, 3,
                                                                       93027, 62502, 94182,
                                                                       41877, 42537, 65802,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 99033, 0, 3,
                                                                       94182, 63327, 95337,
                                                                       42537, 43197, 66792,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 100419, 0, 3,
                                                                       95337, 64152, 96492,
                                                                       43197, 43857, 67782,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101805, 3, 45177,
                                                                       45192, 68814, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101833, 3, 45192,
                                                                       45207, 68835, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101861, 3, 45207,
                                                                       45222, 68856, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101889, 3, 45222,
                                                                       45237, 68877, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101917, 3, 45237,
                                                                       45252, 68898, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101945, 3, 45252,
                                                                       45267, 68919, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101973, 3, 45267,
                                                                       45282, 68940, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102001, 3, 45282,
                                                                       45297, 68961, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102029, 3, 45297,
                                                                       45312, 68982, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102057, 3, 45312,
                                                                       45327, 69003, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102085, 3, 45327,
                                                                       45342, 69024, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102113, 0, 3,
                                                                       101805, 68814, 101833,
                                                                       45372, 45417, 69171,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102197, 0, 3,
                                                                       101833, 68835, 101861,
                                                                       45417, 45462, 69234,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102281, 0, 3,
                                                                       101861, 68856, 101889,
                                                                       45462, 45507, 69297,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102365, 0, 3,
                                                                       101889, 68877, 101917,
                                                                       45507, 45552, 69360,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102449, 0, 3,
                                                                       101917, 68898, 101945,
                                                                       45552, 45597, 69423,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102533, 0, 3,
                                                                       101945, 68919, 101973,
                                                                       45597, 45642, 69486,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102617, 0, 3,
                                                                       101973, 68940, 102001,
                                                                       45642, 45687, 69549,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102701, 0, 3,
                                                                       102001, 68961, 102029,
                                                                       45687, 45732, 69612,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102785, 0, 3,
                                                                       102029, 68982, 102057,
                                                                       45732, 45777, 69675,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102869, 0, 3,
                                                                       102057, 69003, 102085,
                                                                       45777, 45822, 69738,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 102953, 0, 3,
                                                                       102113, 69171, 102197,
                                                                       45912, 46002, 70053,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103121, 0, 3,
                                                                       102197, 69234, 102281,
                                                                       46002, 46092, 70179,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103289, 0, 3,
                                                                       102281, 69297, 102365,
                                                                       46092, 46182, 70305,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103457, 0, 3,
                                                                       102365, 69360, 102449,
                                                                       46182, 46272, 70431,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103625, 0, 3,
                                                                       102449, 69423, 102533,
                                                                       46272, 46362, 70557,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103793, 0, 3,
                                                                       102533, 69486, 102617,
                                                                       46362, 46452, 70683,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103961, 0, 3,
                                                                       102617, 69549, 102701,
                                                                       46452, 46542, 70809,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 104129, 0, 3,
                                                                       102701, 69612, 102785,
                                                                       46542, 46632, 70935,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 104297, 0, 3,
                                                                       102785, 69675, 102869,
                                                                       46632, 46722, 71061,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 104465, 0, 3,
                                                                       102953, 70053, 103121,
                                                                       46902, 47052, 71607,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 104745, 0, 3,
                                                                       103121, 70179, 103289,
                                                                       47052, 47202, 71817,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 105025, 0, 3,
                                                                       103289, 70305, 103457,
                                                                       47202, 47352, 72027,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 105305, 0, 3,
                                                                       103457, 70431, 103625,
                                                                       47352, 47502, 72237,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 105585, 0, 3,
                                                                       103625, 70557, 103793,
                                                                       47502, 47652, 72447,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 105865, 0, 3,
                                                                       103793, 70683, 103961,
                                                                       47652, 47802, 72657,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 106145, 0, 3,
                                                                       103961, 70809, 104129,
                                                                       47802, 47952, 72867,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 106425, 0, 3,
                                                                       104129, 70935, 104297,
                                                                       47952, 48102, 73077,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 106705, 0, 3,
                                                                       104465, 71607, 104745,
                                                                       48402, 48627, 73917,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 107125, 0, 3,
                                                                       104745, 71817, 105025,
                                                                       48627, 48852, 74232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 107545, 0, 3,
                                                                       105025, 72027, 105305,
                                                                       48852, 49077, 74547,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 107965, 0, 3,
                                                                       105305, 72237, 105585,
                                                                       49077, 49302, 74862,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 108385, 0, 3,
                                                                       105585, 72447, 105865,
                                                                       49302, 49527, 75177,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 108805, 0, 3,
                                                                       105865, 72657, 106145,
                                                                       49527, 49752, 75492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 109225, 0, 3,
                                                                       106145, 72867, 106425,
                                                                       49752, 49977, 75807,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 109645, 0, 3,
                                                                       106705, 73917, 107125,
                                                                       50427, 50742, 77004,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 110233, 0, 3,
                                                                       107125, 74232, 107545,
                                                                       50742, 51057, 77445,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 110821, 0, 3,
                                                                       107545, 74547, 107965,
                                                                       51057, 51372, 77886,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 111409, 0, 3,
                                                                       107965, 74862, 108385,
                                                                       51372, 51687, 78327,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 111997, 0, 3,
                                                                       108385, 75177, 108805,
                                                                       51687, 52002, 78768,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 112585, 0, 3,
                                                                       108805, 75492, 109225,
                                                                       52002, 52317, 79209,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 113173, 0, 3,
                                                                       109645, 77004, 110233,
                                                                       52947, 53367, 80826,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 113957, 0, 3,
                                                                       110233, 77445, 110821,
                                                                       53367, 53787, 81414,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 114741, 0, 3,
                                                                       110821, 77886, 111409,
                                                                       53787, 54207, 82002,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 115525, 0, 3,
                                                                       111409, 78327, 111997,
                                                                       54207, 54627, 82590,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 116309, 0, 3,
                                                                       111997, 78768, 112585,
                                                                       54627, 55047, 83178,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 117093, 0, 3,
                                                                       113173, 80826, 113957,
                                                                       55887, 56427, 85278,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 118101, 0, 3,
                                                                       113957, 81414, 114741,
                                                                       56427, 56967, 86034,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 119109, 0, 3,
                                                                       114741, 82002, 115525,
                                                                       56967, 57507, 86790,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 120117, 0, 3,
                                                                       115525, 82590, 116309,
                                                                       57507, 58047, 87546,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 121125, 0, 3,
                                                                       117093, 85278, 118101,
                                                                       59127, 59802, 90192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 122385, 0, 3,
                                                                       118101, 86034, 119109,
                                                                       59802, 60477, 91137,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 123645, 0, 3,
                                                                       119109, 86790, 120117,
                                                                       60477, 61152, 92082,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 124905, 0, 3,
                                                                       121125, 90192, 122385,
                                                                       62502, 63327, 95337,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 126445, 0, 3,
                                                                       122385, 91137, 123645,
                                                                       63327, 64152, 96492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 127985, 0, 3,
                                                                       124905, 95337, 126445,
                                                                       65802, 66792, 100419,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129833, 3, 68772,
                                                                       68793, 101805, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129869, 3, 68793,
                                                                       68814, 101833, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129905, 3, 68814,
                                                                       68835, 101861, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129941, 3, 68835,
                                                                       68856, 101889, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129977, 3, 68856,
                                                                       68877, 101917, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130013, 3, 68877,
                                                                       68898, 101945, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130049, 3, 68898,
                                                                       68919, 101973, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130085, 3, 68919,
                                                                       68940, 102001, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130121, 3, 68940,
                                                                       68961, 102029, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130157, 3, 68961,
                                                                       68982, 102057, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130193, 3, 68982,
                                                                       69003, 102085, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130229, 0, 3,
                                                                       129833, 101805, 129869,
                                                                       69045, 69108, 102113,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130337, 0, 3,
                                                                       129869, 101833, 129905,
                                                                       69108, 69171, 102197,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130445, 0, 3,
                                                                       129905, 101861, 129941,
                                                                       69171, 69234, 102281,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130553, 0, 3,
                                                                       129941, 101889, 129977,
                                                                       69234, 69297, 102365,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130661, 0, 3,
                                                                       129977, 101917, 130013,
                                                                       69297, 69360, 102449,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130769, 0, 3,
                                                                       130013, 101945, 130049,
                                                                       69360, 69423, 102533,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130877, 0, 3,
                                                                       130049, 101973, 130085,
                                                                       69423, 69486, 102617,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130985, 0, 3,
                                                                       130085, 102001, 130121,
                                                                       69486, 69549, 102701,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 131093, 0, 3,
                                                                       130121, 102029, 130157,
                                                                       69549, 69612, 102785,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 131201, 0, 3,
                                                                       130157, 102057, 130193,
                                                                       69612, 69675, 102869,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 131309, 0, 3,
                                                                       130229, 102113, 130337,
                                                                       69801, 69927, 102953,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 131525, 0, 3,
                                                                       130337, 102197, 130445,
                                                                       69927, 70053, 103121,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 131741, 0, 3,
                                                                       130445, 102281, 130553,
                                                                       70053, 70179, 103289,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 131957, 0, 3,
                                                                       130553, 102365, 130661,
                                                                       70179, 70305, 103457,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 132173, 0, 3,
                                                                       130661, 102449, 130769,
                                                                       70305, 70431, 103625,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 132389, 0, 3,
                                                                       130769, 102533, 130877,
                                                                       70431, 70557, 103793,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 132605, 0, 3,
                                                                       130877, 102617, 130985,
                                                                       70557, 70683, 103961,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 132821, 0, 3,
                                                                       130985, 102701, 131093,
                                                                       70683, 70809, 104129,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 133037, 0, 3,
                                                                       131093, 102785, 131201,
                                                                       70809, 70935, 104297,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 133253, 0, 3,
                                                                       131309, 102953, 131525,
                                                                       71187, 71397, 104465,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 133613, 0, 3,
                                                                       131525, 103121, 131741,
                                                                       71397, 71607, 104745,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 133973, 0, 3,
                                                                       131741, 103289, 131957,
                                                                       71607, 71817, 105025,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 134333, 0, 3,
                                                                       131957, 103457, 132173,
                                                                       71817, 72027, 105305,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 134693, 0, 3,
                                                                       132173, 103625, 132389,
                                                                       72027, 72237, 105585,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 135053, 0, 3,
                                                                       132389, 103793, 132605,
                                                                       72237, 72447, 105865,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 135413, 0, 3,
                                                                       132605, 103961, 132821,
                                                                       72447, 72657, 106145,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 135773, 0, 3,
                                                                       132821, 104129, 133037,
                                                                       72657, 72867, 106425,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 136133, 0, 3,
                                                                       133253, 104465, 133613,
                                                                       73287, 73602, 106705,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 136673, 0, 3,
                                                                       133613, 104745, 133973,
                                                                       73602, 73917, 107125,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 137213, 0, 3,
                                                                       133973, 105025, 134333,
                                                                       73917, 74232, 107545,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 137753, 0, 3,
                                                                       134333, 105305, 134693,
                                                                       74232, 74547, 107965,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 138293, 0, 3,
                                                                       134693, 105585, 135053,
                                                                       74547, 74862, 108385,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 138833, 0, 3,
                                                                       135053, 105865, 135413,
                                                                       74862, 75177, 108805,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 139373, 0, 3,
                                                                       135413, 106145, 135773,
                                                                       75177, 75492, 109225,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 139913, 0, 3,
                                                                       136133, 106705, 136673,
                                                                       76122, 76563, 109645,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 140669, 0, 3,
                                                                       136673, 107125, 137213,
                                                                       76563, 77004, 110233,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 141425, 0, 3,
                                                                       137213, 107545, 137753,
                                                                       77004, 77445, 110821,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 142181, 0, 3,
                                                                       137753, 107965, 138293,
                                                                       77445, 77886, 111409,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 142937, 0, 3,
                                                                       138293, 108385, 138833,
                                                                       77886, 78327, 111997,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 143693, 0, 3,
                                                                       138833, 108805, 139373,
                                                                       78327, 78768, 112585,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 144449, 0, 3,
                                                                       139913, 109645, 140669,
                                                                       79650, 80238, 113173,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 145457, 0, 3,
                                                                       140669, 110233, 141425,
                                                                       80238, 80826, 113957,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 146465, 0, 3,
                                                                       141425, 110821, 142181,
                                                                       80826, 81414, 114741,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 147473, 0, 3,
                                                                       142181, 111409, 142937,
                                                                       81414, 82002, 115525,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 148481, 0, 3,
                                                                       142937, 111997, 143693,
                                                                       82002, 82590, 116309,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 149489, 0, 3,
                                                                       144449, 113173, 145457,
                                                                       83766, 84522, 117093,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 150785, 0, 3,
                                                                       145457, 113957, 146465,
                                                                       84522, 85278, 118101,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 152081, 0, 3,
                                                                       146465, 114741, 147473,
                                                                       85278, 86034, 119109,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 153377, 0, 3,
                                                                       147473, 115525, 148481,
                                                                       86034, 86790, 120117,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 154673, 0, 3,
                                                                       149489, 117093, 150785,
                                                                       88302, 89247, 121125,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 156293, 0, 3,
                                                                       150785, 118101, 152081,
                                                                       89247, 90192, 122385,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 157913, 0, 3,
                                                                       152081, 119109, 153377,
                                                                       90192, 91137, 123645,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 159533, 0, 3,
                                                                       154673, 121125, 156293,
                                                                       93027, 94182, 124905,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 161513, 0, 3,
                                                                       156293, 122385, 157913,
                                                                       94182, 95337, 126445,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 163493, 0, 3,
                                                                       159533, 124905, 161513,
                                                                       97647, 99033, 127985,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 165869, 144449, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 167297, 149489, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 169133, 154673, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 171428, 159533, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 174233, 163493, 2376, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 166877, 165869, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 168593, 167297, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 170753, 169133, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 173408, 171428, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 176609, 174233, 66, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 177599, 166877, 168593, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 178859, 168593, 170753, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 180479, 170753, 173408, 15, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 182504, 173408, 176609, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 184979, 177599, 178859, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 187499, 178859, 180479, 15, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 190739, 180479, 182504, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 194789, 184979, 187499, 15, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 198989, 187499, 190739, 15, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 204389, 194789, 198989, 15, nmax);

        simdtrf::transform_i_inner(buffer, 210689, 204389, 15, 15, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 210689, 195, nmax);
    }

    for (size_t m = 0; m < 1755; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
